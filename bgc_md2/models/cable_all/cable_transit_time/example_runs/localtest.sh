a=1
b=2
func(){
  local a b
  a=$((a + 10))
  b=$((a + 10))
  echo $a $b
}
echo $a $b
func $a $b
echo $a $b
