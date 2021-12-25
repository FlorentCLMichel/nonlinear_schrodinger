test_cpu:
	cargo test --offline

test_cpu_multithread:
	cargo test --offline --features "multithread_ft"

test_opencl:
	cargo test --offline --no-default-features --features "opencl"

clippy:
	cargo clippy --offline

debug_cpu:
	cargo build --offline

debug_cpu_multithread:
	cargo build --offline --features "multithread_ft"

debug_opencl:
	cargo build --offline --no-default-features --features "opencl"

release_cpu:
	cargo build --release --offline

release_cpu_multithread:
	cargo build --release --offline --features "multithread_ft"

release_opencl:
	cargo build --release --no-default-features --offline --features "opencl"

clean:
	rm -r target; rm Cargo.lock
