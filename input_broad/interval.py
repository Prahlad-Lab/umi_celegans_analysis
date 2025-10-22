import glob, os
    # Works around a JES limitation where multiples files with the same name overwrite each other when globbed
intervals = sorted(glob.glob("out/*/*.interval_list"))
for i, interval in enumerate(intervals):
    (directory, filename) = os.path.split(interval)
    newName = os.path.join(directory, str(i + 1) + filename)
    os.rename(interval, newName)
print(len(intervals))
if len(intervals) == 0:
    raise ValueError("Interval list produced 0 scattered interval lists. Is the gtf or input interval list empty?")
f = open("interval_count.txt", "w+")
f.write(str(len(intervals)))
f.close()
