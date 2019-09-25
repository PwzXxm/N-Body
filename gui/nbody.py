import pygame
import sys

G_SCREEN_SIZE = [640,480]
G_FPS_LIMIT = 60
G_SHOW_FPS = False



def draw(dt: int, screen: pygame.Surface) -> None:
    pass

def draw_fps(dt: int, screen: pygame.Surface) -> None:
    font = pygame.font.Font(None, 20)
    fps = font.render(str(int(1000 / dt)), True, pygame.Color('white'))
    screen.blit(fps, (50, 50))


def main():
    pygame.init() 
    
    screen = pygame.display.set_mode(G_SCREEN_SIZE) 
    pygame.display.set_caption('N-Body')

    clock = pygame.time.Clock()

    while True:
        dt = clock.tick(G_FPS_LIMIT)

        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                pygame.quit()
                sys.exit()

        screen.fill(pygame.Color('black'))
        draw(dt, screen)
        if G_SHOW_FPS:
            draw_fps(dt, screen)

        pygame.display.update()

if __name__ == "__main__":
    main()
