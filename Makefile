

docker-build:
	docker-compose build

docker-up:
	docker-compose up

docker-exec:
	docker exec -it inv-man /bin/bash

docker-run-main:
	docker exec -t inv-man julia InvestmentManagement/Main.jl
