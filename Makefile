.PHONY: all all_flags clean prog

all: test1.o f.o
	@gcc -m32 -o prog test1.o f.o -lm

test1.o: test1.c
	@gcc -m32 -c -o test1.o test1.c -lm
f.o: f.asm
	@nasm -f elf32 -o f.o f.asm
prog: test1.o f.o
	@gcc -m32 -o prog test1.o f.o -lm


all_flags: test1.o f.o
	@gcc -m32 -o prog test1.o f.o -lm
	@./prog -help -constants -test_root_auto -steps -test_integral_auto
	@rm -rf prog *.o

clean:
	@rm -rf prog *.o
