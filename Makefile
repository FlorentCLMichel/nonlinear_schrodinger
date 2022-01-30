test:
ifdef test
	cargo test $(test) --offline
else 
	cargo test --offline
endif

test_multithread:
ifdef test
	cargo test $(test) --offline --features "multithread_ft"
else 
	cargo test --offline --features "multithread_ft"
endif

test_opencl:
ifdef test
	cargo test $(test) --offline --no-default-features --features "opencl"
else 
	cargo test --offline --features --no-default-features --features "opencl"
endif

clippy:
	cargo clippy --offline

clippy_opencl:
	cargo clippy --offline --features "opencl"

debug:
	cargo build --offline

debug_multithread:
	cargo build --offline --features "multithread_ft"

debug_opencl:
	cargo build --offline --no-default-features --features "opencl"

release:
	cargo build --release --offline

release_multithread:
	cargo build --release --offline --features "multithread_ft"

release_opencl:
	cargo build --release --no-default-features --offline --features "opencl"

clean:
	rm -r target; rm Cargo.lock
