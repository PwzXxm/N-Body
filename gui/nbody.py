import arcade
import collections
import random
import timeit
import sys
import data_reader

G_SCREEN_WIDTH = 1000
G_SCREEN_HEIGHT = 650
G_SCREEN_TITLE = "N-Body Simulation"
G_FPS_LIMIT = 60
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


class MainWindow(arcade.Window):
    """
    Main application class.
    """

    def __init__(self) -> None:

        # Call the parent class and set up the window
        super().__init__(G_SCREEN_WIDTH, G_SCREEN_HEIGHT, G_SCREEN_TITLE)
        self.particle_list = None
        self.set_update_rate(1 / G_FPS_LIMIT)
        self.start_timer = 0
        self.in_reader = None
        self.out_reader = None

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

        self.particle_list = arcade.SpriteList(use_spatial_hash=False)

        for (pos, w) in zip(in_reader.get_init_pos(), in_reader.get_weights()):
            self.particle_list.append(Particle(pos[0], pos[1], w))
        
        p: Particle
        for p in self.particle_list:
            p.change_x = random.uniform(-1, 1)
            p.change_y = random.uniform(-1, 1)

        self.start_timer = timeit.default_timer()


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
        if timeit.default_timer() - self.start_timer > 4:
            p: Particle
            for p in self.particle_list:
                x = p.center_x + p.change_x * delta_time * 40
                y = p.center_y + p.change_y * delta_time * 40
                p.set_position(x, y)

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
            
            self.stat_draw.tick(timeit.default_timer() - start_time)

def main() -> None:
    in_reader = data_reader.InputDataReader(sys.argv[1])
    out_reader = data_reader.OutputDataReader(sys.argv[2])

    if in_reader.n != out_reader.n:
        raise Exception("data files do not match")
    
    window = MainWindow()
    window.setup(in_reader, out_reader)
    arcade.run()

if __name__ == "__main__":
    main()

