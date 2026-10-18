# Human chromosome 19 compilation profile

The human profile uses LLVM 22 and the exact source/compiler fingerprints recorded in `mumento2-human-chr19-llvm22.json`. The profile build helper verifies those fingerprints and the profile digest before compilation.

Training used the first 32 real HPRC chromosome 19 haplotypes from the frozen benchmark collection, k=63, both strands, and one thread. Every one of the 165,121 training matches agrees exactly with independently validated C++ output. Single-thread training avoids contention on the instrumentation runtime's shared counters. An earlier instrumented 512-input/32-thread attempt was stopped for that reason and produced no profile.

The final executable was checked on 205 end-to-end cases and all 512 haplotypes using 32 threads. Timing of the instrumented training run is excluded from the performance comparison. Benchmark JSON records the chosen graph-seed prefix and any process-local allocator preload separately from the compilation profile.

Older bacterial profiles remain available with their original source/compiler fingerprints. The default profile helper now selects this human profile.
