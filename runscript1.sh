
for i in $(seq 0.0 0.1 1.0); do
  for j in $(seq 0.1 0.1 0.4); do
    echo $i $j
    ./shapeIBM 256 1.0 1.0 0.75 0.2 $i 0.5 50.0 10.0 20 $j 1 0.1 3000 100 0 0 1 >> /dev/null 2>&1 & sleep 10;
  done
  wait
  for j in $(seq 0.5 0.1 0.8); do
    echo $i $j
    ./shapeIBM 256 1.0 1.0 0.75 0.2 $i 0.5 50.0 10.0 20 $j 1 0.1 3000 100 0 0 1 >> /dev/null 2>&1 & sleep 10;
  done
  wait
  for j in $(seq 0.9 0.1 1.2); do
    echo $i $j
    ./shapeIBM 256 1.0 1.0 0.75 0.2 $i 0.5 50.0 10.0 20 $j 1 0.1 3000 100 0 0 1 >> /dev/null 2>&1 & sleep 10;
  done
  wait
  #for j in $(seq 1.3 0.1 1.6); do
  #  echo $i $j
  #  ./shapeIBM 256 1.0 1.0 0.75 0.2 $i 0.5 50.0 10.0 20 $j 1 0.1 1500 100 0 0 1 >> /dev/null 2>&1 & sleep 10;
  #done
  #wait
done
