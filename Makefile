all:
	cargo run 

gitaddall:
	git add src 

loc:
	find src -name '*.rs' | xargs wc -l
