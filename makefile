include .env
export $(shell sed 's/=.*//' .env)

SRC := $(wildcard src/*.c)
OBJ := $(subst .c,.o,$(SRC))
ENGINE := engine
CFLAGS := -Wall -O3
CINCLD := -I../openssl/ -I../openssl/crypto/include -I../openssl/include -I/lib -I/engine
CNF_CFLAGS := -pthread -m64 -Wa,--noexecstack
LIB_CFLAGS := -fPIC $(CNF_CFLAGS) $(CFLAGS) -iquote $(_ROOT)/openssl/include/

$(ENGINE)/engineInterface.o : $(ENGINE)/engine.so $(ENGINE)/engineInterface.cpp
	@g++ -o $(ENGINE)/engineInterface.o $(ENGINE)/engineInterface.cpp -L/usr/local/lib/ -lssl -lecvrf -lcrypto
	@engine/engineInterface.o

$(ENGINE)/engine.so : $(ENGINE)/engine.o
	@gcc -shared -o $@ $^ -lssl -lcrypto
	@openssl engine -t -c $(_ROOT)/algo/engine/engine.so

#engine.o : engine.c $(OBJ)
#	@gcc  $(CINCLD) $(LIB_CFLAGS) -o $@ -c $^ -lssl -lcrypto

$(ENGINE)/engine.o : $(ENGINE)/engine.c $(SRC)
	@gcc  $(CINCLD) $(LIB_CFLAGS) -o $@ -c engine/engine.c -lssl -lcrypto

#%.o : %.c
#	@echo what
#	@echo $^
#	@gcc  $(CINCLD) $(LIB_CFLAGS) -o $@ -c $^ -lssl -lcrypto

#%.o : %../.c
#	@gcc  $(CINCLD) $(LIB_CFLAGS) -o $*.o -c $*.c -lssl -lcrypto
#	@gcc -M $*.c > $*.d


clean:
	@rm -f $(ENGINE)/engineInterface.o
	@rm -f $(ENGINE)/engine.o
	@rm -f $(ENGINE)/engine.so

run:
	@$(ENGINE)/engineInterface.o
