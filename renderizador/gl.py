#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

# pylint: disable=invalid-name

"""
Biblioteca Gráfica / Graphics Library.

Desenvolvido por: Sabrina e Matteo
Disciplina: Computação Gráfica
Data:
"""

from math import cos, pi, sin, tan, atan
from PIL.Image import new
import numpy as np
import rotinas
import rotation_matrix
import time         # Para operações com tempo
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
        GL.anim_key = 0
        GL.scale_matrix = np.array([
            [(GL.width/2), 0, 0, 0],
            [0, (GL.height/2), 0, 0],
            [0, 0, 1, 0],
            [0, 0, 0, 1],
        ])

        GL.translator_matrix = np.array([
            [1, 0, 0, 1],
            [0, 1, 0, 1],
            [0, 0, 1, 0],
            [0, 0, 0, 1],
        ])
        
        GL.mirror_matrix = np.array([
            [1, 0, 0, 0],
            [0, -1, 0, 0],
            [0, 0, 1, 0],
            [0, 0, 0, 1],
        ])

        GL.catmull_rom = np.array([
        [-1/2, 3/2, -3/2, 1/2],
        [1, -5/2, 2, -1/2],
        [-1/2, 0, 1/2, 0],
        [0, 1, 0, 0]
        ])

        GL.value_changed = None

    @staticmethod
    def triangleSet(point, colors, lights=None):
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

        # print("jdinog13i1ofvk efldcm")

        tela = np.matmul(GL.scale_matrix, GL.translator_matrix)
        tela = np.matmul(tela, GL.mirror_matrix)
        
        triangles = [np.array([[point[x], point[x + 3], point[x + 6]], [point[x + 1], point[x + 3 + 1], point[x + 6 + 1]], [point[x + 2], point[x + 3 + 2], point[x + 6 + 2]], [1, 1, 1]]) for x in range(0, len(point), 9)]
        
        transformed_triangles = []
        matriz_transformacao = GL.pilha[0]
        for item in GL.pilha[1:]:
            matriz_transformacao = np.matmul(matriz_transformacao, item)

        for triangulo in triangles:

            triangulo = np.matmul(GL.geometry_transformation, triangulo)

            triangulo = np.matmul(GL.viewpoint_lookat, triangulo)
            
            z_s = triangulo[2, :]

            triangulo = np.matmul(GL.viewpoint_projecao, triangulo)

            triangulo = triangulo/triangulo[-1, :]

            triangulo = np.matmul(tela, triangulo)

            triangulo = triangulo[0:3,:]

            triangulo[2, :] = z_s

            triangulo = triangulo.flatten('F')

            for coord in triangulo:
                transformed_triangles.append(coord)
        
        if lights:
            # print("nfejinfkdfknwkefl")
            rotinas.triangleSet3D_lights(transformed_triangles, colors, lights)
        else:
            rotinas.triangleSet2D(transformed_triangles, colors)


    @staticmethod
    def viewpoint(position, orientation, fieldOfView):
        """Função usada para renderizar (na verdade coletar os dados) de Viewpoint."""
        # Na função de viewpoint você receberá a posição, orientação e campo de visão da
        # câmera virtual. Use esses dados para poder calcular e criar a matriz de projeção
        # perspectiva para poder aplicar nos pontos dos objetos geométricos.

        translacao = np.eye(4)
        translacao[:3,3] = -np.array(position)

        rotacao = rotation_matrix.rotation_matrix(orientation).T

        GL.viewpoint_lookat = np.matmul(rotacao, translacao)

        fovy = 2*atan(tan(fieldOfView/2)*(GL.height/(((GL.height**2) + (GL.width**2))**0.5)))
        top = GL.near * tan(fovy)
        right = top * (GL.width/GL.height)

        GL.viewpoint_projecao = np.array([
            [GL.near/right, 0, 0, 0],
            [0, GL.near/top, 0, 0],
            [0, 0, -(GL.far + GL.near)/(GL.far - GL.near), -(2*GL.far*GL.near)/(GL.far - GL.near)],
            [0, 0, -1, 0],
        ])

    @staticmethod
    def transform_in(translation, scale, rotation):
        """Função usada para renderizar (na verdade coletar os dados) de Transform."""
        # A função transform_in será chamada quando se entrar em um nó X3D do tipo Transform
        # do grafo de cena. Os valores passados são a escala em um vetor [x, y, z]
        # indicando a escala em cada direção, a translação [x, y, z] nas respectivas
        # coordenadas e finalmente a rotação por [x, y, z, t] sendo definida pela rotação
        # do objeto ao redor do eixo x, y, z por t radianos, seguindo a regra da mão direita.
        # Quando se entrar em um nó transform se deverá salvar a matriz de transformação dos
        # modelos do mundo em alguma estrutura de pilha.

        # escala
        escala = np.diag([scale[x] if (x < 3) else 1 for x in range(4)])

        translacao = np.eye(4)
        translacao[:3, 3] = translation

        rotacao = rotation_matrix.rotation_matrix(rotation)

        transformation = np.matmul(translacao, rotacao)
        
        GL.geometry_transformation = np.matmul(transformation, escala)
        GL.pilha.append(GL.geometry_transformation)

    @staticmethod
    def transform_out():
        """Função usada para renderizar (na verdade coletar os dados) de Transform."""
        # A função transform_out será chamada quando se sair em um nó X3D do tipo Transform do
        # grafo de cena. Não são passados valores, porém quando se sai de um nó transform se
        # deverá recuperar a matriz de transformação dos modelos do mundo da estrutura de
        # pilha implementada.

        # O print abaixo é só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.
        GL.geometry_transformation = GL.pilha.pop()


    @staticmethod
    def triangleStripSet(point, stripCount, colors, lights=None):
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

        GL.triangleSet(new_point, colors, lights)

    @staticmethod
    def indexedTriangleStripSet(point, index, colors):
        
        # print("\n"*10 + "DEBUG" + "\n"*10)

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

    @staticmethod
    def box(size, colors):
        """Função usada para renderizar Boxes."""
        # A função box é usada para desenhar paralelepípedos na cena. O Box é centrada no
        # (0, 0, 0) no sistema de coordenadas local e alinhado com os eixos de coordenadas
        # locais. O argumento size especifica as extensões da caixa ao longo dos eixos X, Y
        # e Z, respectivamente, e cada valor do tamanho deve ser maior que zero. Para desenha
        # essa caixa você vai provavelmente querer tesselar ela em triângulos, para isso
        # encontre os vértices e defina os triângulos.

        top = []
        bottom = []

        for y in [0,1]:
            for x in [0,1]:
                for z in [0,1]:
                    lista = [size[0]*x, size[1]*y, size[2]*z]

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
        
        GL.triangleStripSet(coords, [len(coords)], colors)


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

        # print("niroafksdnfasdlkjgnzsdk flxcn,gkjzn,rm gdkjfvs.e,nmdfjkwanlsc,d")

        point = []
        new_colors = []

        for i in coordIndex:
            if i != -1:
                for eixo in range(3):
                    point.append(coord[i*3 + eixo])

        if colorPerVertex:
            for i in colorIndex:
                if i != -1:
                    new_colors.append(i)
        else:
            for i in texCoordIndex:
                if i != -1:
                    new_colors.append(i)


        tela = np.matmul(GL.scale_matrix, GL.translator_matrix)
        tela = np.matmul(tela, GL.mirror_matrix)
        
        triangles = [np.array([[point[x], point[x + 3], point[x + 6]], [point[x + 1], point[x + 3 + 1], point[x + 6 + 1]], [point[x + 2], point[x + 3 + 2], point[x + 6 + 2]], [1, 1, 1]]) for x in range(0, len(point), 9)]

        transformed_triangles = []

        for triangulo in triangles:

            triangulo = np.matmul(GL.geometry_transformation, triangulo)

            triangulo = np.matmul(GL.viewpoint_lookat, triangulo)
            
            z_s = triangulo[2, :]

            triangulo = np.matmul(GL.viewpoint_projecao, triangulo)

            triangulo = triangulo/triangulo[-1, :]

            triangulo = np.matmul(tela, triangulo)

            triangulo = triangulo[0:3,:]

            triangulo[2, :] = z_s

            triangulo = triangulo.flatten('F')

            for coordenada in triangulo:
                transformed_triangles.append(coordenada)

        vertices = transformed_triangles

        if texCoordIndex:
            rotinas.triangleSet3D_tex(vertices, new_colors, current_texture, texCoord)
        elif (colorPerVertex and color):
            print(color)
            print(new_colors)
            rotinas.triangleSet3D_color(vertices, color, new_colors)
        else:
            rotinas.triangleSet3D_white(vertices)


    @staticmethod
    def sphere(radius, colors):
        """Função usada para renderizar Esferas."""
        # A função sphere é usada para desenhar esferas na cena. O esfera é centrada no
        # (0, 0, 0) no sistema de coordenadas local. O argumento radius especifica o
        # raio da esfera que está sendo criada. Para desenha essa esfera você vai
        # precisar tesselar ela em triângulos, para isso encontre os vértices e defina
        # os triângulos.

        # O print abaixo é só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.

        radius *= 1

        long_size = 24
        longitudes = [(theta/long_size)*pi for theta in range(0, (2*long_size) + 1, 1)]

        lat_size = 24
        latitudes = [(bizarrao/lat_size)*pi/2 for bizarrao in range(-lat_size, lat_size + 1, 1)]

        direcao = np.array(GL.lighting["direction"] + [1])

        direcao = np.array(direcao)

        direcao = np.reshape(direcao, (-1, 1))

        direcao = np.reshape(direcao, (-1))[:-1]

        for n in range(len(latitudes) - 1):
            
            pontos = []

            for j in range(len(longitudes)):
                for next_n in range(0, 2):
                    theta, bizarrao = [longitudes[j], latitudes[n + next_n]]
                    x = radius*sin(bizarrao)*cos(theta)
                    y = radius*sin(bizarrao)*sin(theta)
                    z = radius*cos(bizarrao)

                    pontos.append(x)
                    pontos.append(y)
                    pontos.append(z)

            pontos += pontos[:3]

            fix = GL.lighting.copy()
            fix["direction"] = direcao
            
            GL.triangleStripSet(pontos, 0, colors, fix)

    @staticmethod
    def navigationInfo(headlight):
        """Características físicas do avatar do visualizador e do modelo de visualização."""
        # O campo do headlight especifica se um navegador deve acender um luz direcional que
        # sempre aponta na direção que o usuário está olhando. Definir este campo como TRUE
        # faz com que o visualizador forneça sempre uma luz do ponto de vista do usuário.
        # A luz headlight deve ser direcional, ter intensidade = 1, cor = (1 1 1),
        # ambientIntensity = 0,0 e direção = (0 0 −1).

        # O print abaixo é só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.
        print("NavigationInfo : headlight = {0}".format(headlight)) # imprime no terminal

    @staticmethod
    def directionalLight(ambientIntensity, color, intensity, direction):
        """Luz direcional ou paralela."""
        # Define uma fonte de luz direcional que ilumina ao longo de raios paralelos
        # em um determinado vetor tridimensional. Possui os campos básicos ambientIntensity,
        # cor, intensidade. O campo de direção especifica o vetor de direção da iluminação
        # que emana da fonte de luz no sistema de coordenadas local. A luz é emitida ao
        # longo de raios paralelos de uma distância infinita.

        # O print abaixo é só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.
        GL.lighting = {
            "ambientIntensity" : ambientIntensity,
            "color" : color,
            "intensity" : intensity,
            "direction" : direction
        }
        # print("DirectionalLight : ambientIntensity = {0}".format(ambientIntensity))
        # print("DirectionalLight : color = {0}".format(color)) # imprime no terminal
        # print("DirectionalLight : intensity = {0}".format(intensity)) # imprime no terminal
        # print("DirectionalLight : direction = {0}".format(direction)) # imprime no terminal

    @staticmethod
    def pointLight(ambientIntensity, color, intensity, location):
        """Luz pontual."""
        # Fonte de luz pontual em um local 3D no sistema de coordenadas local. Uma fonte
        # de luz pontual emite luz igualmente em todas as direções; ou seja, é omnidirecional.
        # Possui os campos básicos ambientIntensity, cor, intensidade. Um nó PointLight ilumina
        # a geometria em um raio de sua localização. O campo do raio deve ser maior ou igual a
        # zero. A iluminação do nó PointLight diminui com a distância especificada.

        # O print abaixo é só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.
        print("PointLight : ambientIntensity = {0}".format(ambientIntensity))
        print("PointLight : color = {0}".format(color)) # imprime no terminal
        print("PointLight : intensity = {0}".format(intensity)) # imprime no terminal
        print("PointLight : location = {0}".format(location)) # imprime no terminal

    @staticmethod
    def fog(visibilityRange, color):
        """Névoa."""
        # O nó Fog fornece uma maneira de simular efeitos atmosféricos combinando objetos
        # com a cor especificada pelo campo de cores com base nas distâncias dos
        # vários objetos ao visualizador. A visibilidadeRange especifica a distância no
        # sistema de coordenadas local na qual os objetos são totalmente obscurecidos
        # pela névoa. Os objetos localizados fora de visibilityRange do visualizador são
        # desenhados com uma cor de cor constante. Objetos muito próximos do visualizador
        # são muito pouco misturados com a cor do nevoeiro.

        # O print abaixo é só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.
        print("Fog : color = {0}".format(color)) # imprime no terminal
        print("Fog : visibilityRange = {0}".format(visibilityRange))

    @staticmethod
    def timeSensor(cycleInterval, loop):
        """Gera eventos conforme o tempo passa."""
        # Os nós TimeSensor podem ser usados para muitas finalidades, incluindo:
        # Condução de simulações e animações contínuas; Controlar atividades periódicas;
        # iniciar eventos de ocorrência única, como um despertador;
        # Se, no final de um ciclo, o valor do loop for FALSE, a execução é encerrada.
        # Por outro lado, se o loop for TRUE no final de um ciclo, um nó dependente do
        # tempo continua a execução no próximo ciclo. O ciclo de um nó TimeSensor dura
        # cycleInterval segundos. O valor de cycleInterval deve ser maior que zero.

        # Deve retornar a fração de tempo passada em fraction_changed

        # O print abaixo é só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.
        print("TimeSensor : cycleInterval = {0}".format(cycleInterval)) # imprime no terminal
        print("TimeSensor : loop = {0}".format(loop))

        # Esse método já está implementado para os alunos como exemplo
        epoch = time.time()  # time in seconds since the epoch as a floating point number.
        fraction_changed = (epoch % cycleInterval) / cycleInterval

        return fraction_changed

    @staticmethod
    def splinePositionInterpolator(set_fraction, key, keyValue, closed):
        """Interpola não linearmente entre uma lista de vetores 3D."""

        if (key[(GL.anim_key + 1)] < set_fraction):
            GL.anim_key = (GL.anim_key + 1)

        if closed:
            if ((GL.anim_key == len(key) - 2) and (key[GL.anim_key] > set_fraction)):
                GL.anim_key = 0

        px = []
        py = []
        pz = []

        for k in range(4):
            idx = (((GL.anim_key + k)%len(key))*3)
            px += [keyValue[idx]]
            py += [keyValue[idx + 1]]
            pz += [keyValue[idx + 2]]

        tzinho = (set_fraction - key[GL.anim_key])/0.2
        t = np.array([(tzinho**p) for p in range(0, 4)])
        t = np.flip(t, 0)
        
        resultado = np.matmul(t, GL.catmull_rom)

        resultadox = np.matmul(resultado, px)
        resultadoy = np.matmul(resultado, py)
        resultadoz = np.matmul(resultado, pz)

        value_changed = [resultadox, resultadoy, resultadoz]

        return value_changed

    @staticmethod
    def orientationInterpolator(set_fraction, key, keyValue):
        """Interpola entre uma lista de valores de rotação especificos."""
        # Interpola rotações são absolutas no espaço do objeto e, portanto, não são cumulativas.
        # Uma orientação representa a posição final de um objeto após a aplicação de uma rotação.
        # Um OrientationInterpolator interpola entre duas orientações calculando o caminho mais
        # curto na esfera unitária entre as duas orientações. A interpolação é linear em
        # comprimento de arco ao longo deste caminho. Os resultados são indefinidos se as duas
        # orientações forem diagonalmente opostas. O campo keyValue possui uma lista com os
        # valores a serem interpolados, key possui uma lista respectiva de chaves
        # dos valores em keyValue, a fração a ser interpolada vem de set_fraction que varia de
        # zeroa a um. O campo keyValue deve conter exatamente tantas rotações 3D quanto os
        # quadros-chave no key.

        # O print abaixo é só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.
        print("OrientationInterpolator : set_fraction = {0}".format(set_fraction))
        print("OrientationInterpolator : key = {0}".format(key)) # imprime no terminal
        print("OrientationInterpolator : keyValue = {0}".format(keyValue))

        # Abaixo está só um exemplo de como os dados podem ser calculados e transferidos
        value_changed = [0, 0, 1, 0]

        return value_changed

    # Para o futuro (Não para versão atual do projeto.)
    def vertex_shader(self, shader):
        """Para no futuro implementar um vertex shader."""

    def fragment_shader(self, shader):
        """Para no futuro implementar um fragment shader."""
