from typing import *
import struct

MAGIC_NUMBER = 9036

class OutputDataReader():
    def __init__(self, output_file: str) -> None:
        self.f = open(output_file, "rb")
        if self._read_int() != MAGIC_NUMBER:
            raise Exception('invalid data file')
        self.n = self._read_int()
        self.m = self._read_int()
        self.cur_m = 0

        
    def __del__(self) -> None:
        self.f.close()

    def _read_int(self) -> int:
        data = self.f.read(4)
        return struct.unpack("i", data)[0]

    def _read_float(self) -> float:
        data = self.f.read(4)
        return struct.unpack("f", data)[0]
    
    def has_data(self):
        return self.cur_m < self.m

    def get_one_step(self) -> List[Tuple[float, float]]:
        if not self.has_data():
            raise Exception("out of data")
        self.cur_m += 1
        return [(self._read_float(), self._read_float()) for _ in range(self.n)]




# reader = OutputDataReader("./../tmp.out")
# print(reader.n)
# print(reader.m)
# print(reader.get_one_step())

