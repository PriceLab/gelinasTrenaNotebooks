build:
	docker build -t gelinas_skin_trena .

bash:
	docker run -ti --rm gelinas_skin_trena bash

shared:
	docker volume create myTest
	docker run -ti --rm -v myTest:/home/trena/sharedData gelinas_skin_trena bash

packages:
	(cd ./pkgs; R CMD build --no-build-vignettes ~/github/BrowserViz-1)
	(cd ./pkgs; R CMD build --no-build-vignettes ~/github/trenaViz)
	(cd ./pkgs; R CMD build --no-build-vignettes ~/github/trena)

