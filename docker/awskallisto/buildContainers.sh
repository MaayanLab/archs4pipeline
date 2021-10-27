# get rid of old stuff
docker rmi -f $(docker images | grep "^<none>" | awk "{print $3}")
docker rm $(docker ps -q -f status=exited)

docker build -f DockerKallisto -t maayanlab/awskallisto .

docker push maayanlab/awskallisto
