embox_funcs.o:embox_funcs.c embox_funcs.h grid.h
	gcc -Wall embox_funcs.c -g -c

initialise.o:initialise.c initialise.h particles.h grid.h embox_funcs.h
	gcc -Wall initialise.c -g -c

embox.o:embox.c particles.h grid.h update.h initialise.h embox_funcs.h
	gcc -Wall embox.c -g -c

update.o:update.c update.h grid.h particles.h
	gcc -Wall update.c -g -c

embox:embox.o update.o initialise.o embox_funcs.o
	gcc embox.o initialise.o embox_funcs.o update.o  -g -lm -o embox


clean:
	-rm embox  *.o
