
all:
	python3 setup.py install --user

clean:
	rm -fr build/
	rm -fr phanotate.egg-info/
	pip3 uninstall -y phanotate
