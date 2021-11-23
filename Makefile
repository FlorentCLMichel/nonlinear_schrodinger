test:
	cargo test --offline --all-features

debug_cpu:
	cargo build --offline

debug_opencl:
	cargo build --offline --features "opencl"

release_cpu:
	cargo build --offline

release_opencl:
	cargo build --offline --features "opencl"

clean:
	rm -r target; rm Cargo.lock
