from typing import *
import arcade
import collections
import random
import timeit
import sys
import data_reader

G_SCREEN_WIDTH = 1000
G_SCREEN_HEIGHT = 650
G_SCREEN_TITLE = "N-Body Simulation"
G_FPS_LIMIT = 24
G_SHOW_STAT = True



class AvgTimeCounter():
    def __init__(self) -> None:
        self.times = collections.deque(maxlen=20)
    
    def tick(self, time: float) -> None:
        self.times.append(time)

    def get_avg(self) -> float:
        total_time = sum(self.times)
        if total_time == 0:
            return 0
        else:
            return total_time / len(self.times)

class Particle(arcade.Sprite):
    def __init__(self, x: float, y: float, scale: float) -> None:
        super().__init__(filename='p.png', center_x=x, center_y=y, scale=scale)
        self.alpha = 200

class ScaleHelper():
    W_MIN_SCALE = 0.2
    W_MAX_SCALE = 0.8
    P_EMPTY_RANGE = 0.1

    def __init__(self, in_reader: data_reader.InputDataReader) -> None:
        # mass
        ws = in_reader.get_masses()
        self.w_min = min(ws)
        self.w_range = max(ws) - min(ws)

        # position
        pos = in_reader.get_init_pos()
        max_x_abs = max(map(lambda x: abs(x[0]), pos))
        max_y_abs = max(map(lambda x: abs(x[1]), pos))
        x_scale = G_SCREEN_WIDTH / (max_x_abs * 2)
        y_scale = G_SCREEN_HEIGHT / (max_y_abs * 2)
        self.pos_scale = min(x_scale, y_scale)
        self.pos_scale *= (1 - self.P_EMPTY_RANGE)

    
    def scale_mass(self, w: float) -> float:
        return (w - self.w_min) * (self.W_MAX_SCALE - self.W_MIN_SCALE) / self.w_range + self.W_MIN_SCALE

    def scale_pos(self, pos: Tuple[float, float]) -> Tuple[float, float]:
        return (pos[0] * self.pos_scale, pos[1] * self.pos_scale)


class MainWindow(arcade.Window):
    """
    Main application class.
    """

    def __init__(self) -> None:

        # Call the parent class and set up the window
        super().__init__(G_SCREEN_WIDTH, G_SCREEN_HEIGHT, G_SCREEN_TITLE)
        self.particle_list = None
        self.set_update_rate(1 / (G_FPS_LIMIT + 1))
        self.start_timer = 0
        self.in_reader = None
        self.out_reader: data_reader.OutputDataReader = None
        self.scale_helper: ScaleHelper = None
        self.playing = False

        if G_SHOW_STAT:
            self.stat_fps = AvgTimeCounter()
            self.stat_update = AvgTimeCounter()
            self.stat_draw = AvgTimeCounter()

    def setup(self, in_reader: data_reader.InputDataReader, out_reader: data_reader.OutputDataReader) -> None:
        """ Set up the game here. Call this function to restart the game. """
        # make (0, 0) the center point
        arcade.set_viewport(-G_SCREEN_WIDTH//2, G_SCREEN_WIDTH//2, -G_SCREEN_HEIGHT//2, G_SCREEN_HEIGHT//2)
        arcade.set_background_color(arcade.csscolor.BLACK)

        self.in_reader = in_reader
        self.out_reader = out_reader
        
        self.scale_helper = ScaleHelper(in_reader)

        self.particle_list = arcade.SpriteList(use_spatial_hash=False)

        for (pos, w) in zip(in_reader.get_init_pos(), in_reader.get_masses()):
            w = self.scale_helper.scale_mass(w)
            pos = self.scale_helper.scale_pos(pos)
            self.particle_list.append(Particle(pos[0], pos[1], w))
        
        p: Particle
        for p in self.particle_list:
            p.change_x = random.uniform(-1, 1)
            p.change_y = random.uniform(-1, 1)

        self.start_timer = timeit.default_timer()
        self.playing = False


    def update(self, delta_time: float) -> None:
        """
        All the logic to move, and the game logic goes here.
        Normally, you'll call update() on the sprite lists that
        need it.
        """
        start_time = None
        if G_SHOW_STAT:
            start_time = timeit.default_timer()
            self.stat_fps.tick(delta_time)
        # update
        if self.start_timer != None and timeit.default_timer() - self.start_timer > 4:
            print("start playing")
            self.playing = True
            self.start_timer = None

        if self.playing:
            if self.out_reader.has_data():
                p: Particle
                for (p, pos) in zip(self.particle_list, self.out_reader.get_one_step()):
                    pos = self.scale_helper.scale_pos(pos)
                    p.set_position(pos[0], pos[1])
            else:
                self.playing = False
                print("stop playing")

        # ------
        if G_SHOW_STAT:
            self.stat_update.tick(timeit.default_timer() - start_time)


        

    def on_draw(self) -> None:
        """ Render the screen. """
        start_time = None
        if G_SHOW_STAT:
            start_time = timeit.default_timer()
        arcade.start_render()
        # draw

        self.particle_list.draw()    

        # ------
        if G_SHOW_STAT:
            f_time = self.stat_fps.get_avg()
            x = 20 - G_SCREEN_WIDTH // 2
            y = G_SCREEN_HEIGHT // 2 - 30
            if f_time != 0:
                arcade.draw_text( "fps: {}".format(int(1 / f_time)), x, y, arcade.color.WHITE, 16)
            y -= 20
            arcade.draw_text("update: {:.3f}".format(self.stat_update.get_avg()), x, y, arcade.color.WHITE, 16)
            y -= 20
            arcade.draw_text("draw: {:.3f}".format(self.stat_draw.get_avg()), x, y, arcade.color.WHITE, 16)
            y -= 20
            arcade.draw_text("frame: {}".format(self.out_reader.cur_m), x, y, arcade.color.WHITE, 16)

            self.stat_draw.tick(timeit.default_timer() - start_time)

def main() -> None:
    in_reader = data_reader.InputDataReader(sys.argv[1])
    out_reader = data_reader.OutputDataReader(sys.argv[2])

    if in_reader.n != out_reader.n:
        raise Exception("data files do not match")
    
    print("Total record:", out_reader.m)

    window = MainWindow()
    window.setup(in_reader, out_reader)
    arcade.run()

if __name__ == "__main__":
    main()

