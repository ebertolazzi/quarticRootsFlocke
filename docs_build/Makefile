INDEX = ./sphinx/_build/html/genindex.html

build:
	make clean
	@cd sphinx; make html; cd ..
	@make install

install:
	@cd sphinx; rm -rf ../../docs/*; cp -a _build/html/. ../../docs; cd ..

clean:
	@rm -rf xml-cpp
	@rm -rf sphinx*/_build
	@cd sphinx; make clean; cd ..
