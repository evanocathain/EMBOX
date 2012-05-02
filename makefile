embox:embox.c particles.h grid.h
	gcc embox.c -g -lm -o embox

clean:
	rm embox 
