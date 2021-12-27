proc split_cmds(str: string): seq[string] =

  var list: seq[string] = @[]
  var str1: string
  var mode: int = 0

  str1 = ""
  
  for s in str:
    if s == '\'' and mode == 0:
      str1 = str1 & s
      mode = 1
    elif s == '\'' and mode == 1:
      str1 = str1 & s
      mode = 0
    elif s == ',' and mode == 0:
      list.add(str1)
      str1 = ""
    else:
      str1 = str1 & s
    
  list.add(str1)


  return list


#var s: string ="poly, conc, print \'lol, lol, lol\', conc"
#echo s
#echo '\''
#echo split_cmds(s)

