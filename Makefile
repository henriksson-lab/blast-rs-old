.PHONY: all bench gitaddall loc

all:
	cargo run

bench:
	cargo build --release
	bash bench/run_benchmark.sh

gitaddall:
	git add src

loc:
	find src -name '*.rs' | xargs wc -l
