#include <iostream>
#include <SFML/Graphics.hpp>
#include <cmath>
using namespace std;
using namespace sf;

struct Object {
    double pi = atan(1) * 4.;
    int n = 10;
    double weight = 100000;
    double r;
    long double x = 0;
    long double y = 0;
    long double z = 0;
    long double vx = 0;
    long double vy = 0;
    long double vz = 0;
    long double vxp;
    long double vyp;
    long double vzp;
    long double R;
    double F;
    double G = 0.0000000001;
    double*** sphere;
    void Init() {
        if (n == 0 || weight == 0)
            throw runtime_error("Error\n");
        r = pow(weight, 1. / 3.);
        double x2;
        double y2;
        double z2;
        double r2;
        double ang1;
        double ang2;
        sphere = new double** [2 * n - 1];
        for (int i = 0; i < 2 * n - 1; i++) {
            sphere[i] = new double* [2 * n + 1];
            for (int j = 0; j < 2 * n + 1; j++)
                sphere[i][j] = new double[3];
        }
        for (int i = n - 1; i < 2 * n - 1; i++) {
            sphere[i][0][0] = 0;
            sphere[i][0][1] = 0;
            sphere[i][0][2] = -r;
            sphere[i][2 * n][0] = 0;
            sphere[i][2 * n][1] = 0;
            sphere[i][2 * n][2] = -r;
            sphere[i][n][0] = 0;
            sphere[i][n][1] = 0;
            sphere[i][n][2] = r;
        }
        ang1 = -pi / 2;
        for (int i = 0; i < n - 1; i++) {
            ang1 += pi / n;
            z2 = sin(ang1) * r;
            r2 = cos(ang1) * r;
            ang2 = 0;
            for (int j = 0; j < 2 * n; j++) {
                if (j == 0) {
                    sphere[i][2 * n][0] = cos(ang2) * r2;
                    sphere[i][2 * n][1] = sin(ang2) * r2;
                    sphere[i][2 * n][2] = z2;
                }
                sphere[i][j][0] = cos(ang2) * r2;
                sphere[i][j][1] = sin(ang2) * r2;
                sphere[i][j][2] = z2;
                sphere[j % n + n - 1][i + 1 + (j / n) * (n - i - 1) * 2][0] = cos(ang2) * r2;
                sphere[j % n + n - 1][i + 1 + (j / n) * (n - i - 1) * 2][1] = sin(ang2) * r2;
                sphere[j % n + n - 1][i + 1 + (j / n) * (n - i - 1) * 2][2] = z2;
                ang2 += pi / n;
            }
        }
    }
    void Physics(long double x, long double y, long double z, double weight, double c) {
        vxp = x - this->x;
        vyp = y - this->y;
        vzp = z - this->z;
        if (vxp == 0 && vyp == 0 && vzp == 0)
            return;
        R = sqrt(vxp * vxp + vyp * vyp + vzp * vzp);
        vxp /= R;
        vyp /= R;
        vzp /= R;
        F = (c * c * G * weight * this->weight * abs(weight)) / (R * R * (abs(this->weight) + abs(weight)));
        vxp *= F;
        vyp *= F;
        vzp *= F;
        vx += vxp;
        vy += vyp;
        vz += vzp;
    }
    void Move() {
        x += vx;
        y += vy;
        z += vz;
    }
};
int main()
{
    int b = 1;
    double c = 1;
    int L = 2;
    Object* obj = new Object[L];
    obj[0].n = 10;
    obj[0].weight = 100000000;
    obj[0].x = 0;
    obj[0].y = 0;
    obj[0].z = 0;
    obj[0].vx = 0;
    obj[0].vy = 0;
    obj[0].vz = 0;
    obj[0].Init();

    obj[1].n = 10;
    obj[1].weight = 100000;
    obj[1].x = 10000;
    obj[1].y = 0;
    obj[1].z = 0;
    obj[1].vx = 0;
    obj[1].vy = sqrt(0.1);
    obj[1].vz = 0;
    obj[1].Init();

    /*obj[2].n = 10;
    obj[2].weight = 10000;
    obj[2].x = 10050;
    obj[2].y = 0;
    obj[2].z = 0;
    obj[2].vx = 0;
    obj[2].vy = sqrt(0.1) + sqrt(0.002);
    obj[2].vz = 0;
    obj[2].Init();*/

    int width = 1200;
    int height = 1000;
    double pi = atan(1) * 4.;
    double pang1 = pi / 2.;
    double pang2 = 0;
    long double px = 0;
    long double py = 0;
    long double pz = 0;
    long double x1;
    long double y1;
    long double z1;
    long double x2;
    long double y2;
    long double z2;
    long double l1;
    long double l2;
    double ang1;
    double ang1p;
    double ang2;
    double ang2p;
    long double coef;
    double x;
    double y;
    double pcos;
    double psin;
    bool permit;
    RenderWindow window(VideoMode(width, height), "SFML");
    RectangleShape win(Vector2f(width, height));
    win.setFillColor(Color::Black);
    while (window.isOpen()) {
        Event event;
        Mouse::setPosition(sf::Vector2i(width / 2, height / 2), window);
        while (window.pollEvent(event)) {
            if (event.type == Event::Closed)
                window.close();
            if (event.type == Event::MouseMoved) {
                pang1 += (width / 2 - event.mouseMove.x) * 0.003;
                pang2 += (height / 2 - event.mouseMove.y) * 0.003;
                if (pang1 < 0)
                    pang1 += pi * 2;
                if (pang1 > pi * 2)
                    pang1 -= pi * 2;
                if (pang2 < -pi / 2)
                    pang2 = -pi / 2;
                if (pang2 > pi / 2)
                    pang2 = pi / 2;
            }
        }
        pcos = cos(pang1);
        psin = sin(pang1);
        if (Keyboard::isKeyPressed(Keyboard::Escape))
            window.close();
        if (Keyboard::isKeyPressed(Keyboard::W)) {
            px += pcos * b;
            py += psin * b;
        }
        if (Keyboard::isKeyPressed(Keyboard::S)) {
            px -= pcos * b;
            py -= psin * b;
        }
        if (Keyboard::isKeyPressed(Keyboard::D)) {
            px += psin * b;
            py -= pcos * b;
        }
        if (Keyboard::isKeyPressed(Keyboard::A)) {
            px -= psin * b;
            py += pcos * b;
        }
        if (Keyboard::isKeyPressed(Keyboard::Space))
            pz += b;
        if (Keyboard::isKeyPressed(Keyboard::LShift))
            pz -= b;
        if (Keyboard::isKeyPressed(Keyboard::Left)) {
            c--;
            if (c == 0)
                c = 1;
            else
                for (int i = 0; i < L; i++) {
                    obj[i].vx *= c / (c + 1);
                    obj[i].vy *= c / (c + 1);
                    obj[i].vz *= c / (c + 1);
                }
        }
        if (Keyboard::isKeyPressed(Keyboard::Right)) {
            c++;
            if (c == 1001)
                c = 1000;
            else
                for (int i = 0; i < L; i++) {
                    obj[i].vx *= c / (c - 1);
                    obj[i].vy *= c / (c - 1);
                    obj[i].vz *= c / (c - 1);
                }
        }
        if (Keyboard::isKeyPressed(Keyboard::Up)) {
            b++;
            if (b == 1001)
                b = 1000;
        }
        if (Keyboard::isKeyPressed(Keyboard::Down)) {
            b--;
            if (b == 0)
                b = 1;
        }
        window.clear();
        window.setMouseCursorVisible(false);
        window.draw(win);
        Vertex line1[] =
        {
            Vertex(Vector2f(width / 2 - 5, height / 2)),
            Vertex(Vector2f(width / 2 + 6, height / 2))
        };
        window.draw(line1, 2, Lines);
        Vertex line2[] =
        {
            Vertex(Vector2f(width / 2, height / 2 - 6)),
            Vertex(Vector2f(width / 2, height / 2 + 5))
        };
        window.draw(line2, 2, Lines);
        for (int i = 0; i < L; i++) {
            for (int j = 0; j < L; j++)
                obj[i].Physics(obj[j].x, obj[j].y, obj[j].z, obj[j].weight, c);
            obj[i].Move();
            for (int j = 0; j < obj[i].n * 2 - 1; j++) {
                permit = false;
                for (int k = 0; k < obj[i].n * 2 + 1; k++) {
                    x1 = obj[i].sphere[j][k][0] + obj[i].x - px;
                    y1 = obj[i].sphere[j][k][1] + obj[i].y - py;
                    z1 = obj[i].sphere[j][k][2] + obj[i].z - pz;
                    l1 = sqrt(x1 * x1 + y1 * y1);
                    ang1 = acos(x1 / l1);
                    if (y1 < 0)
                        ang1 = pi * 2 - ang1;
                    ang1p = pang1 - ang1;
                    if (ang1p > pi)
                        ang1p -= pi * 2;
                    if (ang1p < -pi)
                        ang1p += pi * 2;
                    y2 = l1 * abs(cos(ang1p));
                    ang2 = atan(z1 / y2);
                    if (abs(ang1p) > pi / 2)
                        ang2 = pi - ang2;
                    else if (ang2 < 0)
                        ang2 += pi * 2;
                    if (pang2 >= 0)
                        ang2p = ang2 - pang2;
                    else
                        ang2p = ang2 - pang2 - pi * 2;
                    if (ang2p > pi)
                        ang2p -= pi * 2;
                    if (ang2p < -pi)
                        ang2p += pi * 2;
                    if (abs(ang2p) > pi / 2)
                        continue;
                    l2 = sqrt(y2 * y2 + z1 * z1);
                    coef = l2 * cos(abs(ang2p)) / 350;
                    x2 = l1 * sin(ang1p) / coef;
                    z2 = l2 * sin(ang2p) / coef;
                    if (permit) {
                        Vertex line[] =
                        {
                            Vertex(Vector2f(x, y)),
                            Vertex(Vector2f(width / 2 + x2, height / 2 - z2))
                        };
                        window.draw(line, 2, Lines);
                    }
                    x = width / 2 + x2;
                    y = height / 2 - z2;
                    permit = true;
                }
            }
        }
        window.display();
    }
    return 0;
}
