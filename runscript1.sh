

#for i in $(seq 0.5 0.1 1.0); do
  #for j in $(seq 0.4 0.05 0.6); do
  #  echo $j
  #  #./shapeIBM 256 1.0 1.0 0.75 0.2 $i 0.5 50.0 10.0 20 $j 1 0.1 3000 100 0 0 0 >> /dev/null 2>&1 & sleep 10;
  #  ./shapeIBM 256 1.0 1.0 0.75 0.2 0.4 0.5 50.0 10.0 20 $j 1 0.1 3000 100 0 0 0 >> /dev/null 2>&1 & sleep 10;
  #done
  #wait
   for j in $(seq 0.6 0.05 0.8); do
     echo $j
  #   ./shapeIBM 256 1.0 1.0 0.75 0.2 $i 0.5 50.0 10.0 20 $j 1 0.1 3000 100 0 0 0 >> /dev/null 2>&1 & sleep 10;
      ./shapeIBM 256 1.0 1.0 0.75 0.2 0.4 0.5 50.0 10.0 20 $j 1 0.1 3000 100 0 0 0 >> /dev/null 2>&1 & sleep 10;
  done
  wait
  for j in $(seq 0.8 0.05 1.0); do
    echo $j
    #./shapeIBM 256 1.0 1.0 0.75 0.2 $i 0.5 50.0 10.0 20 $j 1 0.1 3000 100 0 0 0 >> /dev/null 2>&1 & sleep 10;
    ./shapeIBM 256 1.0 1.0 0.75 0.2 0.4 0.5 50.0 10.0 20 $j 1 0.1 3000 100 0 0 0 >> /dev/null 2>&1 & sleep 10;
  done
  wait
  for j in $(seq 1.0 0.05 1.2); do
    echo $j
  #   ./shapeIBM 256 1.0 1.0 0.75 0.2 $i 0.5 50.0 10.0 20 $j 1 0.1 3000 100 0 0 0 >> /dev/null 2>&1 & sleep 10;
    ./shapeIBM 256 1.0 1.0 0.75 0.2 0.4 0.5 50.0 10.0 20 $j 1 0.1 3000 100 0 0 0 >> /dev/null 2>&1 & sleep 10;
  done
  wait
  for j in $(seq 1.2 0.05 1.4); do
    echo $j
  #   ./shapeIBM 256 1.0 1.0 0.75 0.2 $i 0.5 50.0 10.0 20 $j 1 0.1 3000 100 0 0 0 >> /dev/null 2>&1 & sleep 10;
    ./shapeIBM 256 1.0 1.0 0.75 0.2 0.4 0.5 50.0 10.0 20 $j 1 0.1 3000 100 0 0 0 >> /dev/null 2>&1 & sleep 10;
  done
  wait
  for j in $(seq 1.4 0.05 1.6); do
    echo $j
  #  ./shapeIBM 256 1.0 1.0 0.75 0.2 $i 0.5 50.0 10.0 20 $j 1 0.1 3000 100 0 0 0 >> /dev/null 2>&1 & sleep 10;
    ./shapeIBM 256 1.0 1.0 0.75 0.2 0.4 0.5 50.0 10.0 20 $j 1 0.1 3000 100 0 0 0 >> /dev/null 2>&1 & sleep 10;
  done
  wait
  ./shapeIBM 256 1.0 1.0 0.75 0.2 0.4 0.5 50.0 10.0 20 1.6 1 0.1 3000 100 0 0 0 >> /dev/null 2>&1 & sleep 10;

#done
