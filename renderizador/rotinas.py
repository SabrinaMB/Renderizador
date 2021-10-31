#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

"""
Rotinas de operação de nós X3D.

Desenvolvido por:
Disciplina: Computação Gráfica
Data:
"""

import gpu          # Simula os recursos de uma GPU
import numpy as np

#################################################################################
# NÃO USAR MAIS ESSE ARQUIVO. AS ROTINAS DEVEM SER IMPLEMENTADAS AGORA NO gl.GL #
#################################################################################

# web3d.org/documents/specifications/19775-1/V3.0/Part01/components/geometry2D.html#Polypoint2D
def polypoint2D(point, colors):
    """Função usada para renderizar Polypoint2D."""
    for i in range(0, int(len(point)), 2):
        gpu.GPU.set_pixel(int(point[i]), int(point[i + 1]), 255 * colors['emissiveColor'][0],
                          255 * colors['emissiveColor'][1], 255 * colors['emissiveColor'][2])

    # Nessa função você receberá pontos no parâmetro point, esses pontos são uma lista
    # de pontos x, y sempre na ordem. Assim point[0] é o valor da coordenada x do
    # primeiro ponto, point[1] o valor y do primeiro ponto. Já point[2] é a
    # coordenada x do segundo ponto e assim por diante. Assuma a quantidade de pontos
    # pelo tamanho da lista e assuma que sempre vira uma quantidade par de valores.
    # O parâmetro colors é um dicionário com os tipos cores possíveis, para o Polypoint2D
    # você pode assumir o desenho dos pontos com a cor emissiva (emissiveColor).

    # cuidado com as cores, o X3D especifica de (0,1) e o Framebuffer de (0,255)


# web3d.org/documents/specifications/19775-1/V3.0/Part01/components/geometry2D.html#Polyline2D
def polyline2D(lineSegments, colors):
    # """Função usada para renderizar Polyline2D."""
    if len(lineSegments)%4 != 0:
        lineSegments += lineSegments[:2]
        

    for i in range(0, int(len(lineSegments)), 4):
        x_1 = lineSegments[i]
        y_1 = lineSegments[i + 1]
        x_2 = lineSegments[i + 2]
        y_2 = lineSegments[i + 3]
        u = x_1
        v = y_1
        
        if (x_2 - x_1) == 0:
            s = 0
        else:
            s = ((y_2 - y_1) / (x_2 - x_1))

        p = 0
        if s > 1:
            p = 1/s
        elif s < -1:
            p = -1/s
        else:
            p = 1


        if x_2 > x_1:
            while x_2 + p > u:
                gpu.GPU.draw_pixels([int(u), int(v)], gpu.GPU.RGB8, [255/2 * colors['emissiveColor'][0],
                          255/2 * colors['emissiveColor'][1], 255/2 * colors['emissiveColor'][2]])
                if s > 1:
                    u = u + 1 / s
                    v = v + 1
                elif s < -1:
                    u = u - 1 / s
                    v = v - 1
                else:
                    u = u + 1
                    v = v + s
        elif x_2 < x_1:
            while x_2 - p < u:
                gpu.GPU.draw_pixels([int(u), int(v)], gpu.GPU.RGB8, [255/2 * colors['emissiveColor'][0],
                          255/2 * colors['emissiveColor'][1], 255/2 * colors['emissiveColor'][2]])
                if s > 1:
                    u = u - 1 / s
                    v = v - 1
                elif s < -1:
                    u = u + 1 / s
                    v = v + 1
                else:
                    u = u - 1
                    v = v - s
        else:
            if y_2 > y_1:
                while y_2 + p > v:
                    gpu.GPU.draw_pixels([int(u), int(v)], gpu.GPU.RGB8, [255/2 * colors['emissiveColor'][0],
                            255/2 * colors['emissiveColor'][1], 255/2 * colors['emissiveColor'][2]])
                    
                    v = v + 1
            else:
                while y_2 - p < v:
                    gpu.GPU.draw_pixels([int(u), int(v)], gpu.GPU.RGB8, [255/2 * colors['emissiveColor'][0],
                            255/2 * colors['emissiveColor'][1], 255/2 * colors['emissiveColor'][2]])
                    
                    v = v - 1

    # Nessa função você receberá os pontos de uma linha no parâmetro lineSegments, esses
    # pontos são uma lista de pontos x, y sempre na ordem. Assim point[0] é o valor da
    # coordenada x do primeiro ponto, point[1] o valor y do primeiro ponto. Já point[2] é
    # a coordenada x do segundo ponto e assim por diante. Assuma a quantidade de pontos
    # pelo tamanho da lista. A quantidade mínima de pontos são 2 (4 valores), porém a
    # função pode receber mais pontos para desenhar vários segmentos. Assuma que sempre
    # vira uma quantidade par de valores.
    # O parâmetro colors é um dicionário com os tipos cores possíveis, para o Polyline2D
    # você pode assumir o desenho das linhas com a cor emissiva (emissiveColor).

def ambient(lights, surface):
    return lights["ambientIntensity"] * np.array(surface["diffuseColor"]) * 0.2#surface["shininess"]

def dot_3d(vec1, vec2):
    vec1 = np.array(vec1)
    # vec1 = vec1/np.linalg.norm(vec1)
    
    vec2 = np.array(vec2)
    # vec2 = vec2/np.linalg.norm(vec2)

    x = abs(np.dot(vec1, vec2))
    # x = np.dot(vec1, vec2)

    return x

def normal(surface):
    ponto1 = np.array(surface[:3])
    ponto2 = np.array(surface[3:6])
    ponto3 = np.array(surface[6:])
    vec1 = ponto3 - ponto1
    vec2 = ponto2 - ponto1

    result = np.cross(np.array(vec1), np.array(vec2))
    if  np.linalg.norm(result) == 0:
        return 0
    return np.divide(result, np.linalg.norm(result))

def diffuse(lights, surface):
    return lights["intensity"] * np.array(surface["diffuseColor"]) * dot_3d(surface["normal"], lights["direction"])

def specular(lights, surface):
    L_v = np.array([0, 0, 1]) + np.array(lights["direction"]) 
    # L = 
    L_v = (L_v / np.linalg.norm(L_v))
    return lights["intensity"] * np.array(surface["specularColor"]) * (dot_3d(surface["normal"], L_v) ** (surface["shininess"]*128))

def irgb(surface, lights, print_f):
    Oergb = np.array(np.array(surface["emissiveColor"]))
    ambient_i = ambient(lights, surface)
    diffuse_i = diffuse(lights, surface)
    specular_i = specular(lights, surface)
    Ilrgb = np.array(lights["color"])

    result = Oergb + (Ilrgb*(ambient_i + diffuse_i + specular_i))

    if False:#print_f:
        print("-"*200)
        print(f"surface : {surface}")
        print(f"lights : {lights}")
        print(f"Oergb = {Oergb}")
        print(f"ambient_i = {ambient_i}")
        print(f"diffuse_i = {diffuse_i}")
        print(f"specular_i = {specular_i}")
        print(f"(ambient_i + diffuse_i + specular_i) = {(ambient_i + diffuse_i + specular_i)}")
        print(f"Ilrgb = {Ilrgb}")
        print(f"(Ilrgb*(ambient_i + diffuse_i + specular_i) = {(Ilrgb*(ambient_i + diffuse_i + specular_i))}")
        print(f"result = {result}")
        print("-"*200)

    return result
    # return  + sum(np.array(lights["color"])*(ambient(lights, surface) + diffuse(lights, surface) + specular(lights, surface))))

def dot(x_A, y_A, x_B, y_B, x_ponto, y_ponto):
    return (y_B - y_A)*x_ponto - (x_B - x_A)*y_ponto + (x_B - x_A)*y_A - (y_B - y_A)*x_A

def inside(x_A, y_A, x_B, y_B, x_C, y_C, x_ponto, y_ponto):
    d1 = dot(x_A, y_A, x_B, y_B, x_ponto, y_ponto)
    d2 = dot(x_B, y_B, x_C, y_C, x_ponto, y_ponto)
    d3 = dot(x_C, y_C, x_A, y_A, x_ponto, y_ponto)
    if d1 >= 0 and d2 >= 0 and d3 >= 0:
    # if dot(x_A, y_A, x_B, y_B, x_ponto, y_ponto) >= 0 and dot(x_B, y_B, x_C, y_C, x_ponto, y_ponto) >= 0 and dot(x_C, y_C, x_A, y_A, x_ponto, y_ponto) >= 0:
        return 1
    else:
        return 0

def inside_fuq(x_A, y_A, x_B, y_B, x_C, y_C, x_ponto, y_ponto):
    d1 = dot(x_A, y_A, x_B, y_B, x_ponto, y_ponto)
    d2 = dot(x_B, y_B, x_C, y_C, x_ponto, y_ponto)
    d3 = dot(x_C, y_C, x_A, y_A, x_ponto, y_ponto)
    # if d1 >= 0 and d2 >= 0 and d3 >= 0:
    if dot(x_A, y_A, x_B, y_B, x_ponto, y_ponto) >= 0 and dot(x_B, y_B, x_C, y_C, x_ponto, y_ponto) >= 0 and dot(x_C, y_C, x_A, y_A, x_ponto, y_ponto) >= 0:
        return [d1/dot(x_A, y_A, x_B, y_B, x_C, y_C), d2/dot(x_B, y_B, x_C, y_C, x_A, y_A), d3/dot(x_C, y_C, x_A, y_A, x_B, y_B)]
    else:
        return None

# web3d.org/documents/specifications/19775-1/V3.0/Part01/components/geometry2D.html#TriangleSet2D
def triangleSet2D(vertices, colors):
    """Função usada para renderizar TriangleSet2D."""
    for i in range(0, len(vertices), 6):
        lista = vertices[i:i+6]
        x_1, y_1, x_2, y_2, x_3, y_3 = lista
        minX = min([lista[i*2] for i in range(3)])
        maxX = max([lista[i*2] for i in range(3)])
        minY = min([lista[1 + i*2] for i in range(3)])
        maxY = max([lista[1 + i*2] for i in range(3)])

        aliasing = 4

        for x in range(int(minX)*aliasing, int(maxX)*aliasing, aliasing):
            for y in range(int(minY)*aliasing, int(maxY)*aliasing, aliasing):
                sumk = 0
                for kx in range(aliasing):
                    for ky in range(aliasing):
                        sumk += inside(x_1, y_1, x_2, y_2, x_3, y_3, (x + kx) / aliasing, (y + ky) / aliasing)
                
                sumk = sumk*255/aliasing

                if sumk > 0:
                    
                    if sum(colors["emissiveColor"]) > 0:
                        gpu.GPU.draw_pixels([int(x/aliasing), int(y/aliasing)], gpu.GPU.RGB8, [sumk * colors['emissiveColor'][0],
                                          sumk * colors['emissiveColor'][1], sumk * colors['emissiveColor'][2]])
                    elif sum(colors["diffuseColor"]) > 0:
                        gpu.GPU.draw_pixels([int(x / aliasing), int(y / aliasing)], gpu.GPU.RGB8,
                                            [sumk * colors['diffuseColor'][0],
                                             sumk * colors['diffuseColor'][1], sumk * colors['diffuseColor'][2]])
                    else:
                        print("color error")

def triangleSet3D_lights(vertices, material, lights, aliasing=4):
    for k in range(0, len(vertices), 9):

        lista = [vertices[k + i] for i in range(0, 9) if (i + 1) % 3 != 0]
        x_1, y_1, x_2, y_2, x_3, y_3 = lista
        minX = min([lista[i*2] for i in range(3)])
        maxX = max([lista[i*2] for i in range(3)])
        minY = min([lista[1 + i*2] for i in range(3)])
        maxY = max([lista[1 + i*2] for i in range(3)])

        material["normal"] = normal(vertices[k:k + 9])

        aliasing = 1
        print_f = False


        
        for x in range(round(minX), round(maxX) + 1, 1):
            for y in range(round(minY), round(maxY) + 1, 1):
        # for x in range(int(minX)*aliasing, int(maxX)*aliasing, aliasing):
        #     for y in range(int(minY)*aliasing, int(maxY)*aliasing, aliasing):

                sumk = 0

                for kx in range(aliasing):
                    for ky in range(aliasing):

                        sksk = inside_fuq(x_1, y_1, x_2, y_2, x_3, y_3, x + (1/(aliasing*2)) + (kx / aliasing), y + (1/(aliasing*2)) + (ky / aliasing))

                        if sksk is not None:
                            sumk += 1


                if sumk > 0:
                    pixel = irgb(material, lights, print_f)*256
                    
                    sumk = (aliasing**2)/sumk

                    if print_f:
                        print(pixel)
                        print(sumk)
                        print_f = False
                    pixel = np.around(pixel/sumk)
                    gpu.GPU.draw_pixels([x, y], gpu.GPU.RGB8, list(pixel[:3]))

def triangleSet3D_white(vertices, aliasing=4):
    for k in range(0, len(vertices), 9):

        lista = [vertices[k + i] for i in range(0, 9) if (i + 1) % 3 != 0]
        x_1, y_1, x_2, y_2, x_3, y_3 = lista
        minX = min([lista[i*2] for i in range(3)])
        maxX = max([lista[i*2] for i in range(3)])
        minY = min([lista[1 + i*2] for i in range(3)])
        maxY = max([lista[1 + i*2] for i in range(3)])

        aliasing = 8
        
        for x in range(int(minX)*aliasing, int(maxX)*aliasing, aliasing):
            for y in range(int(minY)*aliasing, int(maxY)*aliasing, aliasing):

                sumk = 0

                for kx in range(aliasing):
                    for ky in range(aliasing):

                        sksk = inside_fuq(x_1, y_1, x_2, y_2, x_3, y_3, (x + kx) / aliasing, (y + ky) / aliasing)

                        if sksk is not None:
                            sumk += 1


                if sumk > 0:
                    pixel = np.array([255, 255, 255])
                    sumk = int((aliasing**2)/sumk)
                    pixel //= sumk
                    gpu.GPU.draw_pixels([int(x/aliasing), int(y/aliasing)], gpu.GPU.RGB8, list(pixel[:3]))

def triangleSet3D_color(vertices, color, new_colors, aliasing=4):
    for k in range(0, len(vertices), 9):

        lista = [vertices[k + i] for i in range(0, 9) if (i + 1) % 3 != 0]
        x_1, y_1, x_2, y_2, x_3, y_3 = lista
        minX = min([lista[i*2] for i in range(3)])
        maxX = max([lista[i*2] for i in range(3)])
        minY = min([lista[1 + i*2] for i in range(3)])
        maxY = max([lista[1 + i*2] for i in range(3)])

        z0, z1, z2 = [1/vertices[k + i] for i in range(0, 9) if (i + 1) % 3 == 0]

        rgb0 = [color[new_colors[(k//9)*3 + 0]*3 + bolas] for bolas in range(3)]
        rgb1 = [color[new_colors[(k//9)*3 + 1]*3 + bolas] for bolas in range(3)]
        rgb2 = [color[new_colors[(k//9)*3 + 2]*3 + bolas] for bolas in range(3)]

        for x in range(int(minX)*aliasing, int(maxX)*aliasing, aliasing):
            for y in range(int(minY)*aliasing, int(maxY)*aliasing, aliasing):

                sumk = 0
                rgb = [0, 0, 0]

                for kx in range(aliasing):
                    for ky in range(aliasing):

                        sksk = inside_fuq(x_1, y_1, x_2, y_2, x_3, y_3, (x + kx) / aliasing, (y + ky) / aliasing)

                        if sksk is not None:
                            gamma, alpha, beta = sksk

                            Z_zao = 1/((z0*alpha) + (z1*beta) + (z2*gamma))
                            sumk += 1

                            rgb[0] += Z_zao*((rgb0[0]*z0*alpha) + (rgb1[0]*z1*beta) + (rgb2[0]*z2*gamma))
                            rgb[1] += Z_zao*((rgb0[1]*z0*alpha) + (rgb1[1]*z1*beta) + (rgb2[1]*z2*gamma))
                            rgb[2] += Z_zao*((rgb0[2]*z0*alpha) + (rgb1[2]*z1*beta) + (rgb2[2]*z2*gamma))

                if sum(rgb) > 0:
                    rgb = list(255*np.array(rgb)/(aliasing**2))

                    gpu.GPU.draw_pixels([int(x/aliasing), int(y/aliasing)], gpu.GPU.RGB8, rgb)

def triangleSet3D_tex(vertices, new_colors, current_texture, texCoord, aliasing=4):
    for k in range(0, len(vertices), 9):

        lista = [vertices[k + i] for i in range(0, 9) if (i + 1) % 3 != 0]
        x_1, y_1, x_2, y_2, x_3, y_3 = lista
        minX = min([lista[i*2] for i in range(3)])
        maxX = max([lista[i*2] for i in range(3)])
        minY = min([lista[1 + i*2] for i in range(3)])
        maxY = max([lista[1 + i*2] for i in range(3)])

        z0, z1, z2 = [1/vertices[k + i] for i in range(0, 9) if (i + 1) % 3 == 0]

        image = gpu.GPU.load_texture(current_texture[0])
        rgb0 = [texCoord[new_colors[(k // 9) * 3 + 0] * 2 + bolas] for bolas in range(2)]
        rgb1 = [texCoord[new_colors[(k // 9) * 3 + 1] * 2 + bolas] for bolas in range(2)]
        rgb2 = [texCoord[new_colors[(k // 9) * 3 + 2] * 2 + bolas] for bolas in range(2)]

        for x in range(int(minX)*aliasing, int(maxX)*aliasing, aliasing):
            for y in range(int(minY)*aliasing, int(maxY)*aliasing, aliasing):

                sumk = 0
                rgb = [0, 0]

                for kx in range(aliasing):
                    for ky in range(aliasing):

                        sksk = inside_fuq(x_1, y_1, x_2, y_2, x_3, y_3, (x + kx) / aliasing, (y + ky) / aliasing)

                        if sksk is not None:
                            gamma, alpha, beta = sksk

                            Z_zao = 1/((z0*alpha) + (z1*beta) + (z2*gamma))

                            if sumk == 0:
                                rgb[0] += Z_zao*((rgb0[0]*z0*alpha) + (rgb1[0]*z1*beta) + (rgb2[0]*z2*gamma))
                                rgb[1] += Z_zao*((rgb0[1]*z0*alpha) + (rgb1[1]*z1*beta) + (rgb2[1]*z2*gamma))
                            sumk += 1


                if sumk > 0:
                    pixel = image[image.shape[0] - 1 - int(rgb[1] * image.shape[0]), int(rgb[0] * image.shape[1])]
                    sumk = int((aliasing**2)/sumk)

                    pixel //= sumk
                    gpu.GPU.draw_pixels([int(x/aliasing), int(y/aliasing)], gpu.GPU.RGB8, list(pixel[:3]))
    