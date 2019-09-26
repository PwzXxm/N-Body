import arcade
import collections
import random
import timeit

G_SCREEN_WIDTH = 1000
G_SCREEN_HEIGHT = 650
G_SCREEN_TITLE = "N-Body Simulation"
G_FPS_LIMIT = 60
G_SHOW_STAT = True
G_PARTICLE_COLOR = (0, 255, 187, 200)




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
    def __init__(self) -> None:
        super().__init__(filename='p.png', scale=0.2)
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

        if G_SHOW_STAT:
            self.stat_fps = AvgTimeCounter()
            self.stat_update = AvgTimeCounter()
            self.stat_draw = AvgTimeCounter()

    def setup(self) -> None:
        """ Set up the game here. Call this function to restart the game. """
        arcade.set_background_color(arcade.csscolor.BLACK)

        self.particle_list = arcade.SpriteList(use_spatial_hash=False)

        for _ in range(4000):
            self.particle_list.append(Particle())
        
        p: Particle
        for p in self.particle_list:
            p.set_position(random.randint(0, G_SCREEN_WIDTH), random.randint(0, G_SCREEN_HEIGHT))
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
            if f_time != 0:
                arcade.draw_text( "fps: {}".format(int(1 / f_time)), 20, G_SCREEN_HEIGHT - 30, arcade.color.WHITE, 16)
            arcade.draw_text("update: {:.3f}".format(self.stat_update.get_avg()), 20, G_SCREEN_HEIGHT - 50, arcade.color.WHITE, 16)
            arcade.draw_text("draw: {:.3f}".format(self.stat_draw.get_avg()), 20, G_SCREEN_HEIGHT - 70, arcade.color.WHITE, 16)
            
            self.stat_draw.tick(timeit.default_timer() - start_time)

def main() -> None:
    window = MainWindow()
    window.setup()
    arcade.run()

if __name__ == "__main__":
    main()

