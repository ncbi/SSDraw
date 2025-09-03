format:
	black --line-length 79 .

lint:
	flake8 .

docker-build:
	docker build -t ssdraw .

twine-build:
	python3 -m build

twine-check:
	twine check dist/*

twine-upload-test:
	twine upload --repository testpypi dist/*

twine-upload:
	twine upload dist/*

remove-dist:
	rm -rf dist/*