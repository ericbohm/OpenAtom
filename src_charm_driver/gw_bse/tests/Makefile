include ../config.mk

all:
	make -C ../;cd ../input;./tidy;../build/charmrun +p4 ../build/gw_bse ++local test_configuration2
	make -C ../serial_yale;cd ../serial_yale;./gw.x
	make states

pmatrix:
	$(CHARMC) verify_pmatrix.cpp -o verify_pmatrix
	./verify_pmatrix

states:
	$(CHARMC) verify_states.cpp -o verify_states
	./verify_states
