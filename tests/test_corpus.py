"""Named CLI → real GGCAT → locator → compressed-output correctness tests."""
import compression.zstd as zstd
import fcntl
import gzip
import json
import lzma
import math
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

from oracle import verify_output

ROOT = Path(__file__).resolve().parents[1]
BIN = Path(os.environ.get('MUMENTO2_BIN', ROOT/'target/debug/mumento2')).resolve()
GGCAT = Path(os.environ.get('GGCAT_BIN', ROOT/'ggcat/target/debug/ggcat')).resolve()
CASES = json.loads((ROOT/'tests/corpus/cases.json').read_text())


def header(case, c):
    if case.get('header_padding'): return f'c{c} ' + 'x'*case['header_padding']
    return {'same':'same', 'description':f'c{c} annotation with spaces',
            'unicode':f'c{c}_é_菌'}.get(case.get('header'), f'c{c}')


def encode_records(case, genome):
    parts=[]
    for c, sequence in enumerate(genome):
        name=header(case, c)
        if case.get('fastx')=='fastq':
            quality=case.get('quality','I')
            quality=(quality*((len(sequence)+len(quality)-1)//len(quality)))[:len(sequence)]
            parts.append(f'@{name}\n{sequence}\n+\n{quality}\n')
        else:
            width=case.get('wrap', len(sequence) or 1)
            wrapped='\n'.join(sequence[i:i+width] for i in range(0,len(sequence),width))
            parts.append(f'>{name}\n{wrapped}\n')
    text=''.join(parts)
    if not case.get('final_newline',True): text=text.rstrip('\n')
    if case.get('crlf'): text=text.replace('\n','\r\n')
    return text.encode()


def compress(data, codec):
    if codec=='gzip': return gzip.compress(data, mtime=0)
    if codec=='xz': return lzma.compress(data)
    if codec=='zstd': return zstd.compress(data)
    if codec in ('gzip_members','xz_members','zstd_frames'):
        encoder={'gzip_members':gzip.compress, 'xz_members':lzma.compress, 'zstd_frames':zstd.compress}[codec]
        # A deliberately short first decompressor read must not imply EOF.
        return encoder(data[:9])+encoder(data[9:])
    return data


def write_inputs(directory, case, genomes, prefix):
    paths=[]
    for g,genome in enumerate(genomes):
        codec=case.get('compression','plain')
        if codec=='mixed': codec=('gzip','xz','zstd','plain')[g%4]
        extension={'gzip':'gz','xz':'xz','zstd':'zst','gzip_members':'gz',
                   'xz_members':'xz','zstd_frames':'zst','plain':'plain'}[codec]
        path=directory/f'{prefix}-{g}.{"fq" if case.get("fastx")=="fastq" else "fa"}.{extension}'
        path.write_bytes(compress(encode_records(case,genome),codec))
        paths.append(path)
    return paths


def generation(state): return state/(state/'CURRENT').read_text().strip()


class CorpusTests(unittest.TestCase):
    maxDiff=5000
    def shortDescription(self): return None

    def invoke(self, case, output, inputs, extra=(), initial=True, fail=False):
        command=[str(BIN),'--ggcat',str(GGCAT),'-j',str(case.get('threads',3)),
                 '--chunk-bases',str(case.get('chunk_bases',17)),'-m','0.05',
                 '--contig-coords','-o',str(output)]
        if initial:
            command += ['-k',str(case['k'])]
            if case.get('minimum') is not None: command += ['-N',str(case['minimum'])]
            if case.get('m') is not None: command += ['--locate-minimizer-length',str(case['m'])]
            if case.get('forward_only'): command.append('--forward-only')
        command += list(map(str,extra))+list(map(str,inputs))
        run=subprocess.run(command,capture_output=True,text=True,timeout=90)
        diagnostic=f'{case["name"]}: {command}\n{run.stdout[-3500:]}\n{run.stderr[-4500:]}'
        if fail:
            self.assertNotEqual(run.returncode,0,diagnostic)
            self.assertTrue(run.stderr.strip(),diagnostic)
            self.assertIn(case.get('cli_error_contains',case['error_contains']),run.stderr,diagnostic)
        else: self.assertEqual(run.returncode,0,diagnostic)
        return run

    def verify(self,case,output,genomes,minimum,expected):
        text=zstd.decompress(output.read_bytes()).decode()
        result=verify_output(self,text,genomes,case['k'],minimum,case.get('forward_only',False),expected)
        mapping=Path(str(output)+'.contigs.tsv.zst')
        lines=zstd.decompress(mapping.read_bytes()).decode().splitlines()
        self.assertEqual(lines[0],'genome\tinput\tcontig\tname\tlength')
        expected_map={(g,c):(header(case,c).split()[0],len(sequence))
                      for g,genome in enumerate(genomes) for c,sequence in enumerate(genome)}
        actual={}
        for line in lines[1:]:
            g,path,c,name,length=line.split('\t')
            self.assertTrue(path)
            key=(int(g),int(c))
            self.assertNotIn(key,actual)
            actual[key]=(name,int(length))
        self.assertEqual(actual,expected_map)
        # Mumemto-format lengths sidecar (prefix from output path).
        prefix=str(output)
        for suffix in ('.zst','.tsv','.mums'):
            if prefix.endswith(suffix):
                prefix=prefix[:len(prefix)-len(suffix)]
        lengths=Path(prefix+'.lengths')
        self.assertTrue(lengths.is_file(),lengths)
        length_lines=lengths.read_text().splitlines()
        self.assertEqual(sum(1 for line in length_lines if ' * ' in line),len(genomes))
        return result

    def discovery(self,case,directory):
        inputs=write_inputs(directory,case,case['genomes'],'g')
        output=directory/'result.zst'
        times=directory/'timings.tsv'
        protos=directory/'protos.fa.zst'
        run=self.invoke(case,output,inputs,['--timings',times,'--proto-output',protos])
        self.verify(case,output,case['genomes'],case.get('minimum'),case['expected'])
        self.assertIn(f'{len(case["expected"])} MUMs written',run.stderr)
        phases=[line.split('\t') for line in times.read_text().splitlines()]
        self.assertEqual(phases[0],['phase','wall_seconds'])
        names=[p[0] for p in phases[1:]]
        self.assertEqual(len(names),len(set(names)))
        self.assertTrue({'ggcat_build','pipeline_total','proto_index'} <= set(names))
        self.assertTrue(all(math.isfinite(float(t)) and float(t)>=0 for _,t in phases[1:]))
        proto_text=zstd.decompress(protos.read_bytes()).decode()
        for line in proto_text.splitlines():
            if not line.startswith('>'):
                self.assertGreaterEqual(len(line),case['k'])
                self.assertTrue(set(line.upper()) <= set('ACGT'))

    def update(self,case,directory):
        all_genomes=list(case['genomes'])
        inputs=write_inputs(directory,case,all_genomes,'initial')
        state=directory/'state'
        output=directory/'initial.zst'
        self.invoke(case,output,inputs,['--state',state])
        self.verify(case,output,all_genomes,case.get('minimum'),case['expected'])
        minimum=case.get('minimum')
        first_patterns=json.loads((generation(state)/'state.json').read_text())['patterns']
        for path in inputs: path.unlink()
        for i,step in enumerate(case['updates']):
            old_state=state
            pointer=(old_state/'CURRENT').read_bytes()
            inputs=write_inputs(directory,case,step['genomes'],f'new-{i}')
            extra=['--update',old_state]
            if step.get('fork'):
                state=directory/f'fork-{i}'
                extra += ['--state',state]
            if step.get('minimum') is not None:
                minimum=step['minimum']
                extra += ['-N',str(minimum)]
            output=directory/f'updated-{i}.zst'
            run=self.invoke(case,output,inputs,extra,initial=False)
            all_genomes.extend(step['genomes'])
            keys=self.verify(case,output,all_genomes,minimum,step['expected'])
            info=json.loads((generation(state)/'state.json').read_text())
            self.assertEqual(info['genomes'],len(all_genomes))
            self.assertEqual(info['minimum_genomes'],minimum)
            self.assertEqual((info['k'],info['m'],info['forward_only']),
                             (case['k'],case.get('m',min(case['k'],21)),case.get('forward_only',False)))
            self.assertIn(f'{len(step["genomes"])} new source reads',run.stderr)
            if step.get('fork'): self.assertEqual((old_state/'CURRENT').read_bytes(),pointer)
            if case.get('stable_patterns'): self.assertEqual(info['patterns'],first_patterns)
            for path in inputs: path.unlink()
        fresh=directory/'fresh'
        fresh.mkdir()
        inputs=write_inputs(fresh,case,all_genomes,'g')
        rebuilt=fresh/'rebuilt.zst'
        fresh_case=dict(case,minimum=minimum)
        self.invoke(fresh_case,rebuilt,inputs)
        rebuilt_keys=self.verify(case,rebuilt,all_genomes,minimum,case['updates'][-1]['expected'])
        self.assertEqual(keys,rebuilt_keys,'incremental and fresh output differ')

    def error(self,case,directory):
        error=case['error']
        inputs=write_inputs(directory,case,case['genomes'],'g')
        output=directory/'failure-output.zst'
        sentinel=b'previous output must survive a rejected operation'
        output.write_bytes(sentinel)
        basic={'minimum_zero':['-N','0'],'minimum_exceeds_documents':['-N','3'],
               'threads_zero':['-j','0'],'chunk_zero':['--chunk-bases','0'],
               'minimizer_zero':['--locate-minimizer-length','0'],
               'minimizer_exceeds_k':['--locate-minimizer-length','8'],
               'minimizer_exceeds_word':['-k','40','--locate-minimizer-length','30']}
        # Replace defaults rather than duplicate clap options.
        config=dict(case)
        if error in basic:
            extra=basic[error]
            if error=='threads_zero': config['threads']=0; extra=[]
            if error=='chunk_zero': config['chunk_bases']=0; extra=[]
            if error=='minimizer_exceeds_word': config['k']=40; extra=['--locate-minimizer-length','30']
            self.invoke(config,output,inputs,extra,fail=True)
        elif error=='missing_input': self.invoke(case,output,[directory/'absent.fa'],fail=True)
        elif error=='duplicate_input': self.invoke(case,output,[inputs[0],inputs[0]],fail=True)
        elif error=='missing_state': self.invoke(case,output,inputs,['--update',directory/'absent'],initial=False,fail=True)
        elif error=='missing_output_parent':
            self.invoke(case,directory/'absent'/'out.zst',inputs,fail=True)
        else:
            state=directory/'state'
            initial=directory/'initial.zst'
            self.invoke(case,initial,inputs,['--state',state])
            current=generation(state)
            new_inputs=write_inputs(directory,case,[['ACGTTGCAACCTGGA']],'new')
            extra=['--update',state]
            if error=='existing_initial_state': extra=['--state',state]
            elif error=='update_changed_k': extra += ['-k','9']
            elif error=='update_changed_minimizer': extra += ['--locate-minimizer-length','1']
            elif error=='update_changed_strand': extra += ['--forward-only']
            elif error=='update_duplicate_document': new_inputs=[inputs[0]]
            elif error in ('unsupported_state_version','zero_state_genomes'):
                path=current/'state.json'; info=json.loads(path.read_text())
                info['version' if error=='unsupported_state_version' else 'genomes']=999 if error=='unsupported_state_version' else 0
                path.write_text(json.dumps(info))
            elif error=='missing_current': (state/'CURRENT').unlink()
            elif error=='invalid_current_path': (state/'CURRENT').write_text('../escape')
            elif error=='missing_state_manifest': (current/'state.json').unlink()
            elif error=='invalid_state_json': (current/'state.json').write_text('{broken')
            elif error=='invalid_proto_dna': (current/'protos.fa.zst').write_bytes(zstd.compress(b'>invalid\nACNNNNT\n'))
            elif error=='invalid_packed_magic': (current/'genomes/0.zst').write_bytes(zstd.compress(b'badmagic'))
            elif error=='truncated_context_record': (current/'contexts/0/contexts-0-0.zst').write_bytes(zstd.compress(b'bad-record'))
            elif error=='missing_proto_archive': (current/'protos.fa.zst').unlink()
            elif error=='existing_fork_destination':
                fork=directory/'fork'
                self.invoke(case,directory/'fork-initial.zst',inputs,['--state',fork])
                extra += ['--state',fork]
            lock=None
            if error=='locked_state':
                lock=(state/'.lock').open('r+')
                fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
            before=(state/'CURRENT').read_bytes() if (state/'CURRENT').exists() else None
            try:
                self.invoke(case,output,new_inputs,extra,initial=error=='existing_initial_state',fail=True)
            finally:
                if lock: lock.close()
            after=(state/'CURRENT').read_bytes() if (state/'CURRENT').exists() else None
            self.assertEqual(after,before,'rejected update changed the committed generation')
        self.assertEqual(output.read_bytes(),sentinel,'failure replaced an existing output')


def register(case):
    def test(self):
        with tempfile.TemporaryDirectory(prefix=f'm2-{case["name"][:35]}-') as scratch:
            getattr(self,case['kind'])(case,Path(scratch))
    test.__name__='test_'+case['name']
    test.__doc__=case['description']
    setattr(CorpusTests,test.__name__,test)

for case in CASES: register(case)
assert len(CASES)>=100 and len({c['name'] for c in CASES})==len(CASES)

if __name__=='__main__': unittest.main()
