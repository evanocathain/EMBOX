embox.o:embox.c particles.h grid.h updatefields.h
	gcc -Wall embox.c -g -c

updatefields.o:updatefields.c updatefields.h grid.h
	gcc -Wall updatefields.c -g -c

embox:embox.o updatefields.o
	gcc embox.o updatefields.o -g -lm -o embox

clean:
	-rm embox  *.o
