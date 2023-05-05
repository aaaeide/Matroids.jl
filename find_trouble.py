runs = [[],[]]

with open("out.log") as f:
  run = 0
  while line := f.readline():
    line = line.rstrip()
    if line == "---SKILLETEGN---":
      run = 1
    else:
      runs[run].append(line)

i = 0
while True:
  a = runs[0][i]
  b = runs[1][i]

  if a != b:
    print("PROBLEM FOUND!")
    print(a)
    print(b)
    break

  i += 1