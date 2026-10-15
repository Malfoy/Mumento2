use crate::buckets::readers::binary_reader::ChunkDecoder;
use crate::buckets::writers::compressed_binary_writer::COMPRESSED_BUCKET_MAGIC;
use crate::memory_fs::file::reader::FileReader;
use crate::utils::vec_reader::VecReader;

pub struct CompressedStreamDecoder;

const DEFAULT_BUFFER_SIZE: usize = 1024 * 1024;

impl ChunkDecoder for CompressedStreamDecoder {
    const MAGIC_HEADER: &'static [u8; 16] = COMPRESSED_BUCKET_MAGIC;
    type ReaderType =
        VecReader<zstd::stream::read::Decoder<'static, std::io::BufReader<FileReader>>>;

    fn decode_stream(reader: FileReader, _size: u64) -> Self::ReaderType {
        VecReader::new(
            DEFAULT_BUFFER_SIZE,
            zstd::stream::read::Decoder::new(reader)
                .unwrap()
                .single_frame(),
        )
    }

    fn dispose_stream(stream: Self::ReaderType) -> FileReader {
        stream.into_inner().finish().into_inner()
    }
}
