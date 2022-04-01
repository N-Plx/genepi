SUBDIRS = jetset_lib src

all :
	for dir in $(SUBDIRS); do \
	make -C $$dir; \
	done;	

clean : 
	for dir in $(SUBDIRS); do make -C $$dir clean; done;	
