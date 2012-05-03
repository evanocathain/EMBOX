embox.o:embox.c particles.h grid.h update.h
	gcc -Wall embox.c -g -c

update.o:update.c update.h grid.h
	gcc -Wall update.c -g -c

embox:embox.o update.o
	gcc embox.o update.o -g -lm -o embox

clean:
	-rm embox  *.o
