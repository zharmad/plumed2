awk 'BEGIN{
for(i=0;i<1000;i++){
  print "3"
  print "0 0 0"
  print "Ar 0 0 0"
  print "Ar ",rand(), rand(),rand();
  print "Ar ",rand(), rand(),rand();
}
}'
