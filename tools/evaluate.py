import sys
from typing import *
import utils.data_reader as data_reader

UNTRACKED_THRESHOLD = 100

def distance(p1: Tuple[float, float], p2: Tuple[float, float]) -> float:
    return ((p1[0] - p2[0])**2 + (p1[1] - p2[1])**2)**0.5

def main() -> None:
    out_reader1 = data_reader.OutputDataReader(sys.argv[1])
    out_reader2 = data_reader.OutputDataReader(sys.argv[2])

    pos1 = []
    while out_reader1.has_data():
        pos1 = out_reader1.get_one_step()
    pos2 = []
    while out_reader2.has_data():
        pos2 = out_reader2.get_one_step()

    tracking_count = 0
    total_dif = 0.0
    for (p1, p2) in zip(pos1, pos2):
        d = distance(p1, p2)
        if d <= UNTRACKED_THRESHOLD:
            tracking_count += 1
            total_dif += d
    
    print("tracking rate: {:0.2f}%".format(tracking_count * 100 / len(pos1)))
    if tracking_count > 0:
        print("average differnece: {}".format(total_dif / tracking_count))


if __name__ == "__main__":
    main()
