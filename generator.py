import random

def generate_data(n: int) -> None:
    print(n)
    for _ in range(n):
        px = random.uniform(-100, 100)
        py = random.uniform(-100, 100)
        vx = random.uniform(-10, 10)
        vy = random.uniform(-10, 10)
        w = random.uniform(0.01, 2)
        print(px, py, vx, vy, w)


generate_data(30)
