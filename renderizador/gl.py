#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

"""
Biblioteca Gráfica / Graphics Library.

Desenvolvido por: Sabrina e Matteo
Disciplina: Computação Gráfica
Data:
"""

from math import tan, atan
from PIL.Image import new
import numpy as np
import rotinas
import rotation_matrix
import gpu          # Simula os recursos de uma GPU

class GL:
    """Classe que representa a biblioteca gráfica (Graphics Library)."""

    width = 800   # largura da tela
    height = 600  # altura da tela
    near = 0.01   # plano de corte próximo
    far = 1000    # plano de corte distante
    pilha = []

    # def __init__(self):
    #     self.

    @staticmethod
    def setup(width, height, near=0.01, far=1000):
        """Definir parâmetros para camera de razão de aspecto, plano próximo e distante."""
        GL.width = width
        GL.height = height
        GL.near = near
        GL.far = far
        GL.pilha = []

    @staticmethod
    def triangleSet(point, colors):
        print("triangleSet")
        """Função usada para renderizar TriangleSet."""
        # Nessa função você receberá pontos no parâmetro point, esses pontos são uma lista
        # de pontos x, y, e z sempre na ordem. Assim point[0] é o valor da coordenada x do
        # primeiro ponto, point[1] o valor y do primeiro ponto, point[2] o valor z da
        # coordenada z do primeiro ponto. Já point[3] é a coordenada x do segundo ponto e
        # assim por diante.
        # No TriangleSet os triângulos são informados individualmente, assim os três
        # primeiros pontos definem um triângulo, os três próximos pontos definem um novo
        # triângulo, e assim por diante.
        # O parâmetro colors é um dicionário com os tipos cores possíveis, para o TriangleSet
        # você pode assumir o desenho das linhas com a cor emissiva (emissiveColor).
        scale_matrix = np.array([
            [(GL.width/2), 0, 0, 0],
            [0, (GL.height/2), 0, 0],
            [0, 0, 1, 0],
            [0, 0, 0, 1],
        ])

        translator_matrix = np.array([
            [1, 0, 0, 1],
            [0, 1, 0, 1],
            [0, 0, 1, 0],
            [0, 0, 0, 1],
        ])
        
        mirror_matrix = np.array([
            [1, 0, 0, 0],
            [0, -1, 0, 0],
            [0, 0, 1, 0],
            [0, 0, 0, 1],
        ])

        tela = np.matmul(scale_matrix, translator_matrix)
        tela = np.matmul(tela, mirror_matrix)
        
        triangles = [np.array([[point[x], point[x + 3], point[x + 6]], [point[x + 1], point[x + 3 + 1], point[x + 6 + 1]], [point[x + 2], point[x + 3 + 2], point[x + 6 + 2]], [1, 1, 1]]) for x in range(0, len(point), 9)]
        
        transformed_triangles = []
        matriz_transformacao = GL.pilha[0]
        for item in GL.pilha[1:]:
            matriz_transformacao = np.matmul(matriz_transformacao, item)

        for triangulo in triangles:

            triangulo = np.matmul(matriz_transformacao, triangulo)

            triangulo = np.matmul(GL.viewpoint_lookat, triangulo) 

            triangulo = np.matmul(GL.viewpoint_projecao, triangulo)

            triangulo = triangulo/triangulo[-1, :]

            triangulo = np.matmul(tela, triangulo)[0:2,:].flatten('F')

            # print(triangulo)
            
            for coord in triangulo:
                transformed_triangles.append(coord)
        #
        # if sum(colors['emissiveColor']) == 0:
        #     colors['emissiveColor'] = [0, 0, 1]

            # rotinas.triangleSet2D(triangulo, colors)
        # print(transformed_triangles)
        rotinas.triangleSet2D(transformed_triangles, colors)


    @staticmethod
    def viewpoint(position, orientation, fieldOfView):
        print("viewpoint")
        """Função usada para renderizar (na verdade coletar os dados) de Viewpoint."""
        # Na função de viewpoint você receberá a posição, orientação e campo de visão da
        # câmera virtual. Use esses dados para poder calcular e criar a matriz de projeção
        # perspectiva para poder aplicar nos pontos dos objetos geométricos.

        # position = [-8, 1, 1]
        # orientation = [0, 1, 0, -1.57]

        # translacao
        translacao = np.eye(4)
        translacao[:3,3] = -np.array(position)
        #print(translacao)

        # rotacao
        rotacao = rotation_matrix.rotation_matrix(orientation).T
        #print(rotacao)

        GL.viewpoint_lookat = np.matmul(rotacao, translacao)
        #print(transformation)

        fovy = 2*atan(tan(fieldOfView/2)*(GL.height/(((GL.height**2) + (GL.width**2))**0.5)))
        top = GL.near * tan(fovy)
        right = top * (GL.width/GL.height)

        GL.viewpoint_projecao = np.array([
            [GL.near/right, 0, 0, 0],
            [0, GL.near/top, 0, 0],
            [0, 0, -(GL.far + GL.near)/(GL.far - GL.near), -(2*GL.far*GL.near)/(GL.far - GL.near)],
            [0, 0, -1, 0],
        ])

        # # O print abaixo é só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.
        # print("Viewpoint : ", end='')
        # print("position = {0} ".format(position), end='')
        # print("orientation = {0} ".format(orientation), end='')
        # print("fieldOfView = {0} ".format(fieldOfView))

    @staticmethod
    def transform_in(translation, scale, rotation):
        print("transform_in")
        """Função usada para renderizar (na verdade coletar os dados) de Transform."""
        # A função transform_in será chamada quando se entrar em um nó X3D do tipo Transform
        # do grafo de cena. Os valores passados são a escala em um vetor [x, y, z]
        # indicando a escala em cada direção, a translação [x, y, z] nas respectivas
        # coordenadas e finalmente a rotação por [x, y, z, t] sendo definida pela rotação
        # do objeto ao redor do eixo x, y, z por t radianos, seguindo a regra da mão direita.
        # Quando se entrar em um nó transform se deverá salvar a matriz de transformação dos
        # modelos do mundo em alguma estrutura de pilha.


        # translation = [2, 1, 1]
        # scale = [1, 3, 2]
        # rotation = [1, 0, 0, -1.57]

        # escala
        escala = np.diag([scale[x] if (x < 3) else 1 for x in range(4)])
        #print(escala)

        # translacao
        translacao = np.eye(4)
        translacao[:3, 3] = translation
        #print(translacao)

        # rotacao
        rotacao = rotation_matrix.rotation_matrix(rotation)
        #print(rotacao)

        # translacao x rotacao x escala
        transformation = np.matmul(translacao, rotacao)
        GL.geometry_transformation = np.matmul(transformation, escala)
        # print(f"transformacao: \n{transformation}")
        GL.pilha.append(GL.geometry_transformation)
        # # O print abaixo é só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.
        # print("Transform : ", end='')
        # if translation:
        #     print("translation = {0} ".format(translation), end='') # imprime no terminal
        # if scale:
        #     print("scale = {0} ".format(scale), end='') # imprime no terminal
        # if rotation:
        #     print("rotation = {0} ".format(rotation), end='') # imprime no terminal
        # print("")

    @staticmethod
    def transform_out():
        print("transform_out")
        """Função usada para renderizar (na verdade coletar os dados) de Transform."""
        # A função transform_out será chamada quando se sair em um nó X3D do tipo Transform do
        # grafo de cena. Não são passados valores, porém quando se sai de um nó transform se
        # deverá recuperar a matriz de transformação dos modelos do mundo da estrutura de
        # pilha implementada.

        # O print abaixo é só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.
        GL.geometry_transformation = GL.pilha.pop()


    @staticmethod
    def triangleStripSet(point, stripCount, colors):
        """Função usada para renderizar TriangleStripSet."""
        # A função triangleStripSet é usada para desenhar tiras de triângulos interconectados,
        # você receberá as coordenadas dos pontos no parâmetro point, esses pontos são uma
        # lista de pontos x, y, e z sempre na ordem. Assim point[0] é o valor da coordenada x
        # do primeiro ponto, point[1] o valor y do primeiro ponto, point[2] o valor z da
        # coordenada z do primeiro ponto. Já point[3] é a coordenada x do segundo ponto e assim
        # por diante. No TriangleStripSet a quantidade de vértices a serem usados é informado
        # em uma lista chamada stripCount (perceba que é uma lista). Ligue os vértices na ordem,
        # primeiro triângulo será com os vértices 0, 1 e 2, depois serão os vértices 1, 2 e 3,
        # depois 2, 3 e 4, e assim por diante. Cuidado com a orientação dos vértices, ou seja,
        # todos no sentido horário ou todos no sentido anti-horário, conforme especificado.

        ################### THIS WORKS

        new_point = []

        for x in range(0, len(point) - 9, 3):
            if (x/3)%2 == 0:
                for y in range(9):
                    new_point.append(point[x + y])
            else:
                for k in range(2, -1, -1):
                    for y in range(3):
                        new_point.append(point[x + k*3 + y])

        for y in range(6):
            new_point.append(point[len(point) - 6 + y])
        
        for x in range(3):
            new_point.append(point[len(point) - 9 + x])

        #colors['emissiveColor'] = [0, 1, 0]

        GL.triangleSet(new_point, colors)

        ################### THIS WORKS

        
        

        # # O print abaixo é só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.
        # print("TriangleStripSet : pontos = {0} ".format(point), end='')
        # for i, strip in enumerate(stripCount):
        #     print("strip[{0}] = {1} ".format(i, strip), end='')
        # print("")
        # print("TriangleStripSet : colors = {0}".format(colors)) # imprime no terminal as cores

        # # Exemplo de desenho de um pixel branco na coordenada 10, 10
        # gpu.GPU.draw_pixels([10, 10], gpu.GPU.RGB8, [255, 255, 255])  # altera pixel

    @staticmethod
    def indexedTriangleStripSet(point, index, colors):
        """Função usada para renderizar IndexedTriangleStripSet."""
        # A função indexedTriangleStripSet é usada para desenhar tiras de triângulos
        # interconectados, você receberá as coordenadas dos pontos no parâmetro point, esses
        # pontos são uma lista de pontos x, y, e z sempre na ordem. Assim point[0] é o valor
        # da coordenada x do primeiro ponto, point[1] o valor y do primeiro ponto, point[2]
        # o valor z da coordenada z do primeiro ponto. Já point[3] é a coordenada x do
        # segundo ponto e assim por diante. No IndexedTriangleStripSet uma lista informando
        # como conectar os vértices é informada em index, o valor -1 indica que a lista
        # acabou. A ordem de conexão será de 3 em 3 pulando um índice. Por exemplo: o
        # primeiro triângulo será com os vértices 0, 1 e 2, depois serão os vértices 1, 2 e 3,
        # depois 2, 3 e 4, e assim por diante. Cuidado com a orientação dos vértices, ou seja,
        # todos no sentido horário ou todos no sentido anti-horário, conforme especificado.

        ################### THIS WORKS

        new_point = []

        for x in range(0, len(point) - 9, 3):
            if (x/3)%2 == 0:
                for y in range(9):
                    new_point.append(point[x + y])
            else:
                for k in range(2, -1, -1):
                    for y in range(3):
                        new_point.append(point[x + k*3 + y])

        for y in range(6):
            new_point.append(point[len(point) - 6 + y])
        
        for x in range(3):
            new_point.append(point[len(point) - 9 + x])


        if sum(colors['emissiveColor']) == 0:
            colors['emissiveColor'] = [1, 1, 0]

        GL.triangleSet(new_point, colors)

        ################### THIS WORKS

        # # O print abaixo é só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.
        # print("IndexedTriangleStripSet : pontos = {0}, index = {1}".format(point, index))
        # print("IndexedTriangleStripSet : colors = {0}".format(colors)) # imprime as cores

        # # Exemplo de desenho de um pixel branco na coordenada 10, 10
        # gpu.GPU.draw_pixels([10, 10], gpu.GPU.RGB8, [255, 255, 255])  # altera pixel

    @staticmethod
    def box(size, colors):
        """Função usada para renderizar Boxes."""
        # A função box é usada para desenhar paralelepípedos na cena. O Box é centrada no
        # (0, 0, 0) no sistema de coordenadas local e alinhado com os eixos de coordenadas
        # locais. O argumento size especifica as extensões da caixa ao longo dos eixos X, Y
        # e Z, respectivamente, e cada valor do tamanho deve ser maior que zero. Para desenha
        # essa caixa você vai provavelmente querer tesselar ela em triângulos, para isso
        # encontre os vértices e defina os triângulos.


        # size = [4, 1, 4]

        top = []
        bottom = []

        for y in [0,1]:
            for x in [0,1]:
                for z in [0,1]:
                    lista = [size[0]*x, size[1]*y, size[2]*z]
                    # print(lista)
                    for item in lista:
                        if y == 1:
                            top.append(item)
                        else:
                            bottom.append(item)
        
        for i in range(3):
            top.append(top[i])
            bottom.append(bottom[i])

        GL.triangleStripSet(top, [len(top)], colors)
        GL.triangleStripSet(bottom, [len(bottom)], colors)

        ############### THIS TOOO


        ############### THIS WORKS

        coords = []
        for z in range(2):
            if z == 0:
                for x in range(2):
                    for y in range(2):
                        lista = [size[0]*x, size[1]*y, size[2]*z]
                        # print(lista)
                        for item in lista:
                            coords.append(item)
            elif z==1:
                for x in [1, 0]:
                    for y in range(2):
                        lista = [size[0]*x, size[1]*y, size[2]*z]
                        # print(lista)
                        for item in lista:
                            coords.append(item)
        
        for i in range(3):
            coords.append(coords[i])

        for i in range(0, len(coords), 3):
            print(coords[i:i + 3])
        
        GL.triangleStripSet(coords, [len(coords)], colors)

        ############### THIS WORKS

        # # O print abaixo é só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.
        # print("Box : size = {0}".format(size)) # imprime no terminal pontos
        # print("Box : colors = {0}".format(colors)) # imprime no terminal as cores

        # # Exemplo de desenho de um pixel branco na coordenada 10, 10
        # gpu.GPU.draw_pixels([10, 10], gpu.GPU.RGB8, [255, 255, 255])  # altera pixel

    @staticmethod
    def indexedFaceSet(coord, coordIndex, colorPerVertex, color, colorIndex,
                       texCoord, texCoordIndex, colors, current_texture):
        """Função usada para renderizar IndexedFaceSet."""
        # A função indexedFaceSet é usada para desenhar malhas de triângulos. Ela funciona de
        # forma muito simular a IndexedTriangleStripSet porém com mais recursos.
        # Você receberá as coordenadas dos pontos no parâmetro cord, esses
        # pontos são uma lista de pontos x, y, e z sempre na ordem. Assim coord[0] é o valor
        # da coordenada x do primeiro ponto, coord[1] o valor y do primeiro ponto, coord[2]
        # o valor z da coordenada z do primeiro ponto. Já coord[3] é a coordenada x do
        # segundo ponto e assim por diante. No IndexedFaceSet uma lista de vértices é informada
        # em coordIndex, o valor -1 indica que a lista acabou.
        # A ordem de conexão será de 3 em 3 pulando um índice. Por exemplo: o
        # primeiro triângulo será com os vértices 0, 1 e 2, depois serão os vértices 1, 2 e 3,
        # depois 2, 3 e 4, e assim por diante.
        # Adicionalmente essa implementação do IndexedFace aceita cores por vértices, assim
        # se a flag colorPerVertex estiver habilitada, os vértices também possuirão cores
        # que servem para definir a cor interna dos poligonos, para isso faça um cálculo
        # baricêntrico de que cor deverá ter aquela posição. Da mesma forma se pode definir uma
        # textura para o poligono, para isso, use as coordenadas de textura e depois aplique a
        # cor da textura conforme a posição do mapeamento. Dentro da classe GPU já está
        # implementadado um método para a leitura de imagens.

        point = []
        new_colors = []

        for i in coordIndex:
            if i != -1:
                for eixo in range(3):
                    point.append(coord[i*3 + eixo])

        # print(point)

        if colorPerVertex:
            for i in colorIndex:
                if i != -1:
                    new_colors.append(i)
        else:
            for i in texCoordIndex:
                if i != -1:
                    new_colors.append(i)

        scale_matrix = np.array([
            [(GL.width/2), 0, 0, 0],
            [0, (GL.height/2), 0, 0],
            [0, 0, 1, 0],
            [0, 0, 0, 1],
        ])

        translator_matrix = np.array([
            [1, 0, 0, 1],
            [0, 1, 0, 1],
            [0, 0, 1, 0],
            [0, 0, 0, 1],
        ])
        
        mirror_matrix = np.array([
            [1, 0, 0, 0],
            [0, -1, 0, 0],
            [0, 0, 1, 0],
            [0, 0, 0, 1],
        ])

        tela = np.matmul(scale_matrix, translator_matrix)
        tela = np.matmul(tela, mirror_matrix)
        
        triangles = [np.array([[point[x], point[x + 3], point[x + 6]], [point[x + 1], point[x + 3 + 1], point[x + 6 + 1]], [point[x + 2], point[x + 3 + 2], point[x + 6 + 2]], [1, 1, 1]]) for x in range(0, len(point), 9)]
        # print(triangles)
        transformed_triangles = []

        for triangulo in triangles:

            triangulo = np.matmul(GL.geometry_transformation, triangulo)

            triangulo = np.matmul(GL.viewpoint_lookat, triangulo) 

            triangulo = np.matmul(GL.viewpoint_projecao, triangulo)

            triangulo = triangulo/triangulo[-1, :]

            triangulo = np.matmul(tela, triangulo)

            triangulo = triangulo[0:3,:].flatten('F')

            for coordenada in triangulo:
                transformed_triangles.append(coordenada)

        # if sum(colors['emissiveColor']) == 0:
        #     colors['emissiveColor'] = [0, 0, 1]

        vertices = transformed_triangles

        for k in range(0, len(vertices), 9):

            lista = [vertices[k + i] for i in range(0, 9) if (i + 1) % 3 != 0]
            x_1, y_1, x_2, y_2, x_3, y_3 = lista
            minX = min([lista[i*2] for i in range(3)])
            maxX = max([lista[i*2] for i in range(3)])
            minY = min([lista[1 + i*2] for i in range(3)])
            maxY = max([lista[1 + i*2] for i in range(3)])

            # print(lista)

            aliasing = 4
            rgb_print = []

            z0, z1, z2 = [1/vertices[k + i] for i in range(0, 9) if (i + 1) % 3 == 0]

            if colorPerVertex:
                rgb0 = [color[new_colors[(k//9)*3 + 2]*3 + bolas] for bolas in range(3)]
                rgb1 = [color[new_colors[(k//9)*3 + 0]*3 + bolas] for bolas in range(3)]
                rgb2 = [color[new_colors[(k//9)*3 + 1]*3 + bolas] for bolas in range(3)]

                for x in range(int(minX)*aliasing, int(maxX)*aliasing, aliasing):
                    for y in range(int(minY)*aliasing, int(maxY)*aliasing, aliasing):

                        sumk = 0
                        rgb = [0, 0, 0]

                        for kx in range(aliasing):
                            for ky in range(aliasing):

                                sksk = rotinas.inside_fuq(x_1, y_1, x_2, y_2, x_3, y_3, (x + kx) / aliasing, (y + ky) / aliasing)

                                if sksk is not None:
                                    alpha, beta, gamma = sksk

                                    Z_zao = 1/((z0*alpha) + (z1*beta) + (z2*gamma))
                                    sumk += 1

                                    rgb[0] += Z_zao*((rgb0[0]*z0*alpha) + (rgb1[0]*z1*beta) + (rgb2[0]*z2*gamma))
                                    rgb[1] += Z_zao*((rgb0[1]*z0*alpha) + (rgb1[1]*z1*beta) + (rgb2[1]*z2*gamma))
                                    rgb[2] += Z_zao*((rgb0[2]*z0*alpha) + (rgb1[2]*z1*beta) + (rgb2[2]*z2*gamma))

                        # rgb[0] *= 255/aliasing
                        # rgb[1] *= 255/aliasing
                        # rgb[2] *= 255/aliasing

                        if sum(rgb) > 0:
                            rgb = list(255*np.array(rgb)/(aliasing**2))

                            gpu.GPU.draw_pixels([int(x/aliasing), int(y/aliasing)], gpu.GPU.RGB8, rgb)

            else:
                image = gpu.GPU.load_texture(current_texture[0])
                rgb0 = [texCoord[new_colors[(k // 9) * 3 + 0] * 2 + bolas] for bolas in range(2)]
                rgb1 = [texCoord[new_colors[(k // 9) * 3 + 2] * 2 + bolas] for bolas in range(2)]
                rgb2 = [texCoord[new_colors[(k // 9) * 3 + 1] * 2 + bolas] for bolas in range(2)]
                print("texCoord: ", texCoord)
                print("new_colors: ", new_colors)
                print("rgb0: ", rgb0)
                print("rgb1: ", rgb1)
                print("rgb2: ", rgb2)
                for x in range(int(minX)*aliasing, int(maxX)*aliasing, aliasing):
                    for y in range(int(minY)*aliasing, int(maxY)*aliasing, aliasing):

                        sumk = 0
                        rgb = [0, 0]

                        for kx in range(aliasing):
                            for ky in range(aliasing):

                                sksk = rotinas.inside_fuq(x_1, y_1, x_2, y_2, x_3, y_3, (x + kx) / aliasing, (y + ky) / aliasing)

                                if sksk is not None:
                                    alpha, beta, gamma = sksk

                                    Z_zao = 1/((z0*alpha) + (z1*beta) + (z2*gamma))

                                    if sumk == 0:
                                        rgb[0] += Z_zao*((rgb0[0]*z0*alpha) + (rgb1[0]*z1*beta) + (rgb2[0]*z2*gamma))
                                        rgb[1] += Z_zao*((rgb0[1]*z0*alpha) + (rgb1[1]*z1*beta) + (rgb2[1]*z2*gamma))
                                    sumk += 1
                                    # rgb[2] += Z_zao*((rgb0[2]*z0*alpha) + (rgb1[2]*z1*beta) + (rgb2[2]*z2*gamma))
                        #
                        # rgb[0] /= aliasing
                        # rgb[1] /= aliasing
                        # rgb[2] *= 255/aliasing

                        if sumk > 0:
                            pixel = image[int(rgb[0] * image.shape[0]), image.shape[1] - int(rgb[1] * image.shape[1])-1]
                            sumk = int((aliasing**2)/sumk)
                            # print(sumk)
                            pixel //= sumk
                            gpu.GPU.draw_pixels([int(x/aliasing), int(y/aliasing)], gpu.GPU.RGB8, pixel[:3])
        # Os prints abaixo são só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.
        print("IndexedFaceSet : ")
        if coord:
            print("\tpontos(x, y, z) = {0}, coordIndex = {1}".format(coord, coordIndex))
        print("colorPerVertex = {0}".format(colorPerVertex))
        if colorPerVertex and color and colorIndex:
            print("\tcores(r, g, b) = {0}, colorIndex = {1}".format(color, colorIndex))
        if texCoord and texCoordIndex:
            print("\tpontos(u, v) = {0}, texCoordIndex = {1}".format(texCoord, texCoordIndex))
        if current_texture:
            image = gpu.GPU.load_texture(current_texture[0])
            print("\t Matriz com image = {0}".format(image))
            print("\t Dimensões da image = {0}".format(image.shape))
        print("IndexedFaceSet : colors = {0}".format(colors))  # imprime no terminal as cores

        # Exemplo de desenho de um pixel branco na coordenada 10, 10
        gpu.GPU.draw_pixels([10, 10], gpu.GPU.RGB8, [255, 255, 255])  # altera pixel

    # Para o futuro (Não para versão atual do projeto.)
    def vertex_shader(self, shader):
        """Para no futuro implementar um vertex shader."""

    def fragment_shader(self, shader):
        """Para no futuro implementar um fragment shader."""
