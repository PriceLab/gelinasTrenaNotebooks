build:
	docker build -t gelinas_skin_nbserver .

bash:
	docker run -ti --rm gelinas_skin_nbserver bash

shared:
	docker volume create gelinasData
	docker run -ti --rm -v gelinasData:/home/jovyan/sharedData gelinas_skin_nbserver bash

