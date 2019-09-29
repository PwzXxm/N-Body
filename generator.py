import random

def generate_data(n: int) -> None:
    # Unit:
    #  10^5 kg
    #  10^10 m
    print(n)
    print(0.66742)
    for _ in range(n):
        px = random.uniform(-100, 100)
        py = random.uniform(-100, 100)
        vx = random.uniform(-5, 5)
        vy = random.uniform(-5, 5)
        # vx = 0
        # vy = 0
        w = random.uniform(0.01, 20)
        print(px, py, vx, vy, w)


generate_data(100)
