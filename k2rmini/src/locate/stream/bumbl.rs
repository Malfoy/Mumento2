//! Mumemto-compatible `.bumbl` binary MUM writer.
//!
//! Layout (little-endian), matching `mumsio.hpp`:
//! - `uint16` flags: bit15=length32 (always), bit13=partial if any start=-1
//! - `uint64` n_seqs, `uint64` n_mums
//! - `uint32` lengths[n_mums]
//! - `int64` starts[n_mums][n_seqs] (`-1` = absent)
//! - packed strands, MSB-first within each byte
//! - no collinear block section
use super::*;
use std::io::Write;

const FLAG_PARTIAL: u16 = 1 << 13;
const FLAG_LENGTH32: u16 = 1 << 15;

/// One MUM row for bumbl: length, per-genome global starts (`None` → -1), strands.
pub(super) struct BumblRow {
    pub len: u32,
    pub starts: Vec<Option<i64>>,
    pub reverse: Vec<bool>,
}

pub(super) fn write_bumbl(path: &Path, n_seqs: u64, rows: &[BumblRow]) -> io::Result<()> {
    let n_mums = rows.len() as u64;
    let partial = rows
        .iter()
        .any(|r| r.starts.iter().any(Option::is_none));
    let flags = FLAG_LENGTH32 | if partial { FLAG_PARTIAL } else { 0 };
    let parent = path
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or(Path::new("."));
    let mut tmp = tempfile::NamedTempFile::new_in(parent)?;
    {
        let mut out = BufWriter::new(tmp.as_file_mut());
        out.write_all(&flags.to_le_bytes())?;
        out.write_all(&n_seqs.to_le_bytes())?;
        out.write_all(&n_mums.to_le_bytes())?;
        for row in rows {
            out.write_all(&row.len.to_le_bytes())?;
        }
        for row in rows {
            assert_eq!(row.starts.len() as u64, n_seqs);
            for start in &row.starts {
                let v = start.unwrap_or(-1);
                out.write_all(&v.to_le_bytes())?;
            }
        }
        let total = (n_mums as usize).saturating_mul(n_seqs as usize);
        if total > 0 {
            let mut packed = vec![0u8; total.div_ceil(8)];
            let mut idx = 0usize;
            for row in rows {
                assert_eq!(row.reverse.len() as u64, n_seqs);
                for &rev in &row.reverse {
                    // Mumemto: strand bit 1 = forward (+), 0 = reverse (−); MSB-first.
                    if !rev {
                        packed[idx / 8] |= 1 << (7 - (idx % 8));
                    }
                    idx += 1;
                }
            }
            out.write_all(&packed)?;
        }
        out.flush()?;
    }
    tmp.persist(path).map_err(|e| e.error)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn read_u16(b: &[u8], i: &mut usize) -> u16 {
        let v = u16::from_le_bytes(b[*i..*i + 2].try_into().unwrap());
        *i += 2;
        v
    }
    fn read_u64(b: &[u8], i: &mut usize) -> u64 {
        let v = u64::from_le_bytes(b[*i..*i + 8].try_into().unwrap());
        *i += 8;
        v
    }
    fn read_u32(b: &[u8], i: &mut usize) -> u32 {
        let v = u32::from_le_bytes(b[*i..*i + 4].try_into().unwrap());
        *i += 4;
        v
    }
    fn read_i64(b: &[u8], i: &mut usize) -> i64 {
        let v = i64::from_le_bytes(b[*i..*i + 8].try_into().unwrap());
        *i += 8;
        v
    }

    #[test]
    fn bumbl_round_trip_partial_and_strands() -> io::Result<()> {
        let dir = tempfile::tempdir()?;
        let path = dir.path().join("t.bumbl");
        let rows = vec![
            BumblRow {
                len: 40,
                starts: vec![Some(10), Some(20), None],
                reverse: vec![false, true, false],
            },
            BumblRow {
                len: 21,
                starts: vec![Some(0), Some(5), Some(7)],
                reverse: vec![false, false, true],
            },
        ];
        write_bumbl(&path, 3, &rows)?;
        let bytes = std::fs::read(&path)?;
        let mut i = 0;
        let flags = read_u16(&bytes, &mut i);
        assert_eq!(flags & FLAG_LENGTH32, FLAG_LENGTH32);
        assert_eq!(flags & FLAG_PARTIAL, FLAG_PARTIAL);
        assert_eq!(read_u64(&bytes, &mut i), 3);
        assert_eq!(read_u64(&bytes, &mut i), 2);
        assert_eq!(read_u32(&bytes, &mut i), 40);
        assert_eq!(read_u32(&bytes, &mut i), 21);
        assert_eq!(read_i64(&bytes, &mut i), 10);
        assert_eq!(read_i64(&bytes, &mut i), 20);
        assert_eq!(read_i64(&bytes, &mut i), -1);
        assert_eq!(read_i64(&bytes, &mut i), 0);
        assert_eq!(read_i64(&bytes, &mut i), 5);
        assert_eq!(read_i64(&bytes, &mut i), 7);
        let packed = &bytes[i..];
        // 6 bits: + - + (absent) | + + -  → forward bits 1,0,1, 1,1,0
        assert_eq!(packed.len(), 1);
        assert_eq!(packed[0], 0b1011_1000);
        Ok(())
    }
}
