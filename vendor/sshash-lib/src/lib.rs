// SSHash: Sparse and Skew Hashing of k-mers
//
// A Rust implementation of the SSHash compressed k-mer dictionary,
// providing efficient storage and querying of k-mer sets.

#![doc = include_str!("../README.md")]
#![warn(missing_docs)]

pub mod builder;
pub mod constants;
pub(crate) mod darray;
pub mod dictionary;
pub mod encoding;
pub mod hasher;
pub mod kmer;
pub mod kmer_dict_trait;
pub mod minimizer;
pub mod minimizers_control_map;
pub mod mphf_config;
pub mod offsets;
pub mod partitioned_mphf;
pub mod serialization;
pub mod sparse_and_skew_index;
pub mod spectrum_preserving_string_set;
pub mod streaming_query;

/// Packed DNA types used by the SIMD minimizer and query interfaces.
pub use simd_minimizers::packed_seq;
/// Hashing interfaces for the same packed DNA representation.
pub use simd_minimizers::seq_hash;

// Re-export common types at crate root
pub use builder::{BuildConfiguration, CfSegData, DictionaryBuilder, parse_cf_seg};
pub use dictionary::Dictionary;
pub use kmer::{Kmer, Kmer21, Kmer31, Kmer63, KmerBits};
pub use kmer_dict_trait::{KmerDictionary, KmerStreamingQuery};
pub use minimizer::{
    MinimizerInfo, MinimizerIterator, canonical_mmer_at, compute_minimizer, is_self_rc_mmer,
    reverse_complement_mmer,
};
pub use minimizers_control_map::{BucketType, MinimizersControlMap, MinimizersControlMapBuilder};
pub use partitioned_mphf::PartitionedMphf;
pub use streaming_query::{LookupResult, StreamingQuery};

/// Dispatch to the correct const generic `K` based on a runtime `k` value.
///
/// The const generic `K` determines both the storage type (`u64` for K ≤ 31,
/// `u128` for K > 31) **and** the k-mer length used in operations like
/// `reverse_complement()`, `from_str()`, and `decode_kmer()`.  Therefore `K`
/// must always equal the actual k value, not merely the maximum for the
/// storage class.
///
/// `k` must be an odd value in \[3, 63\]; these are the only values for which
/// [`KmerBits`] is implemented.
///
/// # Usage
///
/// ```ignore
/// use sshash_lib::dispatch_on_k;
///
/// dispatch_on_k!(k, K => {
///     let kmer = Kmer::<K>::from_string(s)?;
///     dict.lookup::<K>(&kmer)
/// })
/// ```
///
/// # Panics
///
/// Panics at runtime if `k` is even or outside the \[3, 63\] range.
#[macro_export]
macro_rules! dispatch_on_k {
    ($k:expr, $K:ident => $body:expr) => {{
        match $k {
            3 => {
                const $K: usize = 3;
                $body
            }
            4 => {
                const $K: usize = 4;
                $body
            }
            5 => {
                const $K: usize = 5;
                $body
            }
            6 => {
                const $K: usize = 6;
                $body
            }
            7 => {
                const $K: usize = 7;
                $body
            }
            8 => {
                const $K: usize = 8;
                $body
            }
            9 => {
                const $K: usize = 9;
                $body
            }
            10 => {
                const $K: usize = 10;
                $body
            }
            11 => {
                const $K: usize = 11;
                $body
            }
            12 => {
                const $K: usize = 12;
                $body
            }
            13 => {
                const $K: usize = 13;
                $body
            }
            14 => {
                const $K: usize = 14;
                $body
            }
            15 => {
                const $K: usize = 15;
                $body
            }
            16 => {
                const $K: usize = 16;
                $body
            }
            17 => {
                const $K: usize = 17;
                $body
            }
            18 => {
                const $K: usize = 18;
                $body
            }
            19 => {
                const $K: usize = 19;
                $body
            }
            20 => {
                const $K: usize = 20;
                $body
            }
            21 => {
                const $K: usize = 21;
                $body
            }
            22 => {
                const $K: usize = 22;
                $body
            }
            23 => {
                const $K: usize = 23;
                $body
            }
            24 => {
                const $K: usize = 24;
                $body
            }
            25 => {
                const $K: usize = 25;
                $body
            }
            26 => {
                const $K: usize = 26;
                $body
            }
            27 => {
                const $K: usize = 27;
                $body
            }
            28 => {
                const $K: usize = 28;
                $body
            }
            29 => {
                const $K: usize = 29;
                $body
            }
            30 => {
                const $K: usize = 30;
                $body
            }
            31 => {
                const $K: usize = 31;
                $body
            }
            32 => {
                const $K: usize = 32;
                $body
            }
            33 => {
                const $K: usize = 33;
                $body
            }
            34 => {
                const $K: usize = 34;
                $body
            }
            35 => {
                const $K: usize = 35;
                $body
            }
            36 => {
                const $K: usize = 36;
                $body
            }
            37 => {
                const $K: usize = 37;
                $body
            }
            38 => {
                const $K: usize = 38;
                $body
            }
            39 => {
                const $K: usize = 39;
                $body
            }
            40 => {
                const $K: usize = 40;
                $body
            }
            41 => {
                const $K: usize = 41;
                $body
            }
            42 => {
                const $K: usize = 42;
                $body
            }
            43 => {
                const $K: usize = 43;
                $body
            }
            44 => {
                const $K: usize = 44;
                $body
            }
            45 => {
                const $K: usize = 45;
                $body
            }
            46 => {
                const $K: usize = 46;
                $body
            }
            47 => {
                const $K: usize = 47;
                $body
            }
            48 => {
                const $K: usize = 48;
                $body
            }
            49 => {
                const $K: usize = 49;
                $body
            }
            50 => {
                const $K: usize = 50;
                $body
            }
            51 => {
                const $K: usize = 51;
                $body
            }
            52 => {
                const $K: usize = 52;
                $body
            }
            53 => {
                const $K: usize = 53;
                $body
            }
            54 => {
                const $K: usize = 54;
                $body
            }
            55 => {
                const $K: usize = 55;
                $body
            }
            56 => {
                const $K: usize = 56;
                $body
            }
            57 => {
                const $K: usize = 57;
                $body
            }
            58 => {
                const $K: usize = 58;
                $body
            }
            59 => {
                const $K: usize = 59;
                $body
            }
            60 => {
                const $K: usize = 60;
                $body
            }
            61 => {
                const $K: usize = 61;
                $body
            }
            62 => {
                const $K: usize = 62;
                $body
            }
            63 => {
                const $K: usize = 63;
                $body
            }
            64 => {
                const $K: usize = 64;
                $body
            }
            other => panic!("Unsupported k value: {}. Expected 3..=64.", other),
        }
    }};
}

/// Version information
pub fn version() -> (u8, u8, u8) {
    constants::VERSION
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_version() {
        let (major, minor, patch) = version();
        assert_eq!(major, 0);
        assert_eq!(minor, 1);
        assert_eq!(patch, 0);
    }
}
