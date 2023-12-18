#include "rasterizer.h"
#include "lodepng.h"

#include <Eigen/Dense>
#include <optional>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>

using std::string, Eigen::VectorXd, Eigen::Vector4d, Eigen::Matrix4d;

/* Image data globals */
std::vector<unsigned char> image_;                                                      /* PNG file */
int w_ = -1;                                                                            /* PNG width */
int h_ = -1;                                                                            /* PNG height*/
string png_ = "";                                                                       /* PNG filename */

/* Rasterizer globals */
std::vector<VectorXd> positions_, colors_;                                              /* Buffers */
bool SWAP_ = false;                                                                     /* Swap flag for DDA */

/* Rasterizer elective globals */
std::vector<unsigned int> elements_;                                                    /* Elements elective */
bool depth_ = false;                                                                    /* Depth elective */
std::vector<std::vector<double>> dbuffer_;                                              /* Depth buffer */
Matrix4d mat_;                                                                          /* Uniform matrix */
bool sRGB_ = false;                                                                     /* sRGB elective */
bool hyp_ = false;                                                                      /* Hyperbolic elective */
bool unif_ = false;                                                                     /* Uniform matrix elective*/

std::pair<std::pair<std::vector<unsigned char>, string>, std::pair<int, int>> generateImage(std::ifstream& file) {
    setImageGlobals(file);
    processCommands(file);
    return std::make_pair(std::make_pair(image_, png_), std::make_pair(w_, h_));
}

void drawArraysTriangles(unsigned first, unsigned count) {
    if (count % 3 != 0 ) { std::cerr << "Count must be a multiple of 3\n"; exit(1); }
    for (int i = 0; i < count; i += 3) {
        // std::cout << "drawing triangle with vertices " << first+i << ", " << first+i+1 << ", " << first+i+2 << "\n";
        VectorXd v1(8), v2(8), v3(8);

        // consolidate position & color into 8D vectors
        v1.segment(0,4) = viewportAndDivision(positions_[first+i]), v1.segment(4,4) = divideRGB(colors_[first+i], positions_[first+i][3]);
        v2.segment(0,4) = viewportAndDivision(positions_[first+i+1]), v2.segment(4,4) = divideRGB(colors_[first+i+1], positions_[first+i+1][3]);
        v3.segment(0,4) = viewportAndDivision(positions_[first+i+2]), v3.segment(4,4) = divideRGB(colors_[first+i+2], positions_[first+i+2][3]);

        // std::cout << "input vectors after viewport and w division:\n";
        // std::cout << v1.transpose() << std::endl;
        // std::cout << v2.transpose() << std::endl;
        // std::cout << v3.transpose() << std::endl;

        scanLine(v1, v2, v3);
        SWAP_ = false;
    }
}

void drawElementsTriangles(unsigned count, unsigned offset) {
    if (count % 3 != 0 ) { std::cerr << "Count must be a multiple of 3\n"; exit(1); }
    for (int i = 0; i < count; i += 3) {
        // std::cout << "drawing triangle with elements " << offset+i << ", " << offset+i+1 << ", " << offset+i+2 << "\n";
        VectorXd v1(8), v2(8), v3(8);

        // consolidate position & color into 8D vectors
        v1.segment(0,4) = viewportAndDivision(positions_[elements_[offset+i]]), v1.segment(4,4) = divideRGB(colors_[elements_[offset+i]], positions_[elements_[offset+i]][3]);
        v2.segment(0,4) = viewportAndDivision(positions_[elements_[offset+i+1]]), v2.segment(4,4) = divideRGB(colors_[elements_[offset+i+1]], positions_[elements_[offset+i+1]][3]);
        v3.segment(0,4) = viewportAndDivision(positions_[elements_[offset+i+2]]), v3.segment(4,4) = divideRGB(colors_[elements_[offset+i+2]], positions_[elements_[offset+i+2]][3]);

        // std::cout << "input vectors after viewport and w division:\n";
        // std::cout << v1.transpose() << std::endl;
        // std::cout << v2.transpose() << std::endl;
        // std::cout << v3.transpose() << std::endl;

        scanLine(v1, v2, v3);
        SWAP_ = false;
    }
}

VectorXd divideRGB(VectorXd vec, double w) {
    if (hyp_) {
        // divide RGBA by w
        vec[0] /= w;
        vec[1] /= w;
        vec[2] /= w;
    }
    return vec;
}

VectorXd viewportAndDivision(VectorXd vec) {
    if (unif_) {
        Vector4d result = mat_ * vec.segment(0, 4);
        vec[0] = result[0];
        vec[1] = result[1]; 
        vec[2] = result[2];
        vec[3] = result[3];
    }

    // divide by w
    double w = vec[3];
    vec[0] /= w;
    vec[1] /= w;
    vec[2] /= w;

    if (hyp_) {
        // replace w with 1/w
        vec[3] = (1.0 / w);
    }

    // map depth to [0, 1] range
    if (!hyp_) {
        double depthValue = (vec[2] + 1.0) * 0.5;
        depthValue = std::max(0.0, std::min(depthValue, 1.0));
        vec[3] = depthValue;
    }

    // viewport
    vec[0] = (vec[0] + 1) * w_ / 2;
    vec[1] = (vec[1] + 1) * h_ / 2;
    return vec;
}

void draw(VectorXd &v) {
    if (depth_) {
        // calculate the depth for each pixel
        double pixelDepth = hyp_ ? 1 / v[3] : v[3];
        unsigned int x = v[0];
        unsigned int y = v[1];
        if (pixelDepth < dbuffer_[y][x]) {
            dbuffer_[y][x] = pixelDepth;
        } else {
            return;
        }
    }

    int x = v[0], y = v[1];
    double r = v[4], g = v[5], b = v[6], a = v[7];
    if (hyp_) {
        double w_prime = v[3];
        r /= w_prime; 
        g /= w_prime;
        b /= w_prime;
    }

    // alpha blending (elective)
    double destR = image_[((y * w_) + x) * 4 + 0] / 255.0;
    double destG = image_[((y * w_) + x) * 4 + 1] / 255.0;
    double destB = image_[((y * w_) + x) * 4 + 2] / 255.0;
    double destA = image_[((y * w_) + x) * 4 + 3] / 255.0;
    double blendedR = r + (destR * (1.0 - a));
    double blendedG = g + (destG * (1.0 - a));
    double blendedB = b + (destB * (1.0 - a));
    double blendedA = a + destA * (1.0 - a);
    r = blendedR;
    g = blendedG;
    b = blendedB;
    a = blendedA;

    if (sRGB_) {
        // clamp and convert pixel from RGB to sRGB
        auto clamp = [](double n) -> double { return std::max(0.0, std::min(1.0, n)); };
        auto sRGB = [](double L) -> double { return L > 0.0031308 ? 1.055*pow(L, 1/2.4) - 0.055 : 12.92*L; };
        r = sRGB(clamp(r));
        g = sRGB(clamp(g));
        b = sRGB(clamp(b));
    }

    // draw the pixel
    image_[((y * w_) + x) * 4 + 0] = r * 255.0;
    image_[((y * w_) + x) * 4 + 1] = g * 255.0;
    image_[((y * w_) + x) * 4 + 2] = b * 255.0;
    image_[((y * w_) + x) * 4 + 3] = a * 255.0;

}

std::optional<std::pair<VectorXd, VectorXd>> dda(VectorXd a, VectorXd b, bool axis) {
    unsigned int d = axis ? 1 : 0;
    if (a(d) == b(d)) { 
        return std::nullopt;    // no points found
    }
    if (a(d) > b(d)) {
        if (d == 0) {
            SWAP_ = true;
        } else { 
            std::cout << "undefined behavior: swap on y\n"; 
        }
        VectorXd temp = a;
        a = b;
        b = temp;
    }
    VectorXd delta = b - a;
    VectorXd s = delta / delta(d);
    double e = ceil(a(d)) - a(d);
    VectorXd o = e * s;
    VectorXd p = a + o;
    return std::make_pair(p, s);
}

void scanLine(const VectorXd &P, const VectorXd &Q, const VectorXd &R) {
    // find t, m, b
    VectorXd t = P;
    if (Q[1] < t[1]) {
        t = Q;
    }
    if (R[1] < t[1]) {
        t = R;
    }

    VectorXd b = P;
    if (Q[1] > b[1]) {
        b = Q;
    }
    if (R[1] > b[1]) {
        b = R;
    }

    VectorXd m;
    if (P != t && P != b) {
        m = P;
    } else if (Q != t && Q != b) {
        m = Q;
    } else {
        m = R;
    }

    // std::cout << "t = " << t(1) << " m = " << m(1) << " b = " << b(1) << std::endl;

    // t_b -> top-bot -> longest
    // t_m -> top-mid -> upper half
    // m_b -> mid-bot -> lower half

    auto t_to_b = dda(t, b, 1);
    if (!t_to_b.has_value()) { 
        return;
    }
    VectorXd t_b = t_to_b.value().first;
    VectorXd t_b_step = t_to_b.value().second;

    /* Render top of triangle */
    auto t_to_m = dda(t, m, 1);
    if (t_to_m.has_value()) {
        VectorXd t_m = t_to_m.value().first;
        VectorXd t_m_step = t_to_m.value().second;
        while (t_m(1) < m(1)) {
            double y = t_m(1);
            // std::cout << "y=" << y << "\t";
            // std::cout << "x from " << t_b(0) << " (rgb=" << t_b(4) << ") to " << t_m(0) << " (rgb=" << t_m(4) << ")\t";
            auto tb_to_tm = dda(t_b, t_m, 0);
            if (tb_to_tm.has_value()) {
                VectorXd tb_tm = tb_to_tm.value().first;
                VectorXd tb_tm_step = tb_to_tm.value().second;
                // std::cout << '{';

                double limit1;
                if (SWAP_) {
                    limit1 = t_b(0);
                } else {
                    limit1 = t_m(0);
                }
                SWAP_ = false;

                while (tb_tm(0) < limit1) {
                    double x = tb_tm(0);
                    // std::cout << x << " ";
                    if (y < 0.0 || x < 0.0 || y > static_cast<double>(h_) || x > static_cast<double>(w_)) {
                        // do nothing
                    } else {
                        draw(tb_tm);
                    }
                    tb_tm += tb_tm_step;    // move horizontally (x)
                }
                // std::cout << "}\n";
            }
            t_m += t_m_step;    // move vertically (y)
            t_b += t_b_step;
        }
    }

    /* Render bottom of triangle */
    auto m_to_b = dda(m, b, 1);
    if (m_to_b.has_value()) {
        VectorXd m_b = m_to_b.value().first;
        VectorXd m_b_step = m_to_b.value().second;
        while (m_b(1) < b(1)) {
            double y = m_b(1);
            // std::cout << "y=" << y << "\t";
            // std::cout << "x from " << t_b(0) << " (rgb=" << t_b(4) << ") to " << m_b(0) << " (rgb=" << m_b(4) << ")\t";
            auto tb_to_mb = dda(t_b, m_b, 0);
            if (tb_to_mb.has_value()) {
                VectorXd tb_mb = tb_to_mb.value().first;
                VectorXd tb_mb_step = tb_to_mb.value().second;
                // std::cout << '{';

                double limit2;
                if (SWAP_) {
                    limit2 = t_b(0);
                } else {
                    limit2 = m_b(0);
                }
                SWAP_ = false;

                while (tb_mb(0) < limit2) {
                    double x = tb_mb(0);
                    // std::cout << x << " ";
                    if (y < 0.0 || x < 0.0 || y > static_cast<double>(h_) || x > static_cast<double>(w_)) {
                        // do nothing
                    } else {
                        draw(tb_mb);
                    }
                    tb_mb += tb_mb_step;    // move horizontally (x)
                }
                // std::cout << "}\n";
            }
            m_b += m_b_step;    // move vertically (y)
            t_b += t_b_step;
        }
    }

    return;
}

void processCommands(std::ifstream& file, bool debug) {
    string line;
    int line_n = 1;
    while (std::getline(file, line)) {
        if (!line.empty()) {
            std::istringstream iss(line);
            string cmd;
            iss >> cmd;

            if (cmd == "drawArraysTriangles") {
                string first, count;
                iss >> first >> count;
                drawArraysTriangles(stoi(first), stoi(count));
            } else if (cmd == "drawElementsTriangles") {
                string count, offset;
                iss >> count >> offset;
                drawElementsTriangles(stoi(count), stoi(offset));
            } else if (cmd == "elements") {
                while (!iss.eof()) {
                    string i;
                    iss >> i;
                    elements_.emplace_back(stoi(i));
                }
                std::cout << "elements buffer size " << elements_.size() << std::endl;
            } else if (cmd == "position") {
                string s;
                iss >> s;
                unsigned size = stoi(s);

                std::vector<VectorXd> positions;
                while(!iss.eof()) {
                    VectorXd position(4);
                    switch (size) {
                        case 2: {
                            string x, y;
                            iss >> x >> y;
                            position << stod(x), stod(y), 0.0, 1.0;
                            break;
                        }
                        case 3: {
                            string x, y, z;
                            iss >> x >> y >> z;
                            position << stod(x), stod(y), stod(z), 1.0;
                            break;
                        }
                        case 4: {
                            string x, y, z, w;
                            iss >> x >> y >> z >> w;
                            position << stod(x), stod(y), stod(z), stod(w);
                            break;
                        }
                        default: {
                            std::cerr << "Invalid position size\n";
                            exit(1);
                        }
                    }
                    positions.emplace_back(position);
                }
                positions_ = positions;
            } else if (cmd == "color") {
                string s;
                iss >> s;
                unsigned size = stoi(s);

                std::vector<VectorXd> colors;
                while(!iss.eof()) {
                    VectorXd color(4);
                    switch (size) {
                        case 3: {
                            string r, g, b;
                            iss >> r >> g >> b;
                            color << stod(r), stod(g), stod(b), 1.0;
                            break;
                        }
                        case 4: {
                            string r, g, b, a;
                            iss >> r >> g >> b >> a;
                            color << stod(r), stod(g), stod(b), stod(a);
                            break;
                        }
                        default: {
                            std::cerr << "Invalid color size\n";
                            exit(1);
                        }
                    }
                    colors.emplace_back(color);
                }
                colors_ = colors;
            } else if (cmd == "depth") {
                dbuffer_ = std::vector<std::vector<double>>(h_, std::vector<double>(w_, std::numeric_limits<double>::max()));
                depth_ = true;
            } else if (cmd == "sRGB") {
                sRGB_ = true;
            } else if (cmd == "hyp") {
                hyp_ = true;
            } else if (cmd == "uniformMatrix") {
                string n0,  n1,  n2,  n3,
                       n4,  n5,  n6,  n7,
                       n8,  n9,  n10, n11,
                       n12, n13, n14, n15;
                iss >> n0 >> n1 >> n2 >> n3 >> n4 >> n5 >> n6 >> n7 >> n8 >> n9 >> n10 >> n11 >> n12 >> n13 >> n14 >> n15;
                Matrix4d mat;
                mat(0, 0) = stod(n0);  mat(0, 1) = stod(n4);  mat(0, 2) = stod(n8);   mat(0, 3) = stod(n12);
                mat(1, 0) = stod(n1);  mat(1, 1) = stod(n5);  mat(1, 2) = stod(n9);   mat(1, 3) = stod(n13);
                mat(2, 0) = stod(n2);  mat(2, 1) = stod(n6);  mat(2, 2) = stod(n10);  mat(2, 3) = stod(n14);
                mat(3, 0) = stod(n3);  mat(3, 1) = stod(n7);  mat(3, 2) = stod(n11);  mat(3, 3) = stod(n15);
                mat_ = mat;
                unif_ = true;
                std::cout << mat_ << std::endl;
            }
        }
        line_n++;
    }

    if (debug) {
        std::cout << "positions:\n";
        for (const VectorXd &p : positions_) { std::cout << p.transpose() << std::endl; }
        std::cout << "colors:\n";
        for (const VectorXd &c : colors_) { std::cout << c.transpose() << std::endl; }
    }
}

void setImageGlobals(std::ifstream& file, bool debug) {
    string line;
    std::getline(file, line);

    std::istringstream iss(line);
    string word;
    int i = 0;
    while (iss >> word) {
        if (i != 0) {
            if (i == 1) {
                w_ = stoi(word);
            } else if (i == 2) {
                h_ = stoi(word);
            } else if (i == 3) {
                png_ = word;
            } else {
                std::cerr << "First line had undefined behavior.\n";
                exit(1);
            }
        }
        i++;
    }
    if (debug) { std::cout << "File: " << png_ << " " << w_ << "x" << h_ << '\n'; }
    std::vector<unsigned char> image(w_ * h_ * 4, 0);
    image_ = image;
}

std::ifstream openFile(int argc, char** argv, bool debug) {
    if (argc < 2) {
        std::cerr << "No filename was provided.\n";
        exit(1);
    }
    std::ifstream file(argv[1]);
    if (!file.is_open()) {
        std::cerr << "Could not open file.\n";
        exit(1);
    }
    if (debug) {
        string line;
        std::cout << "File '" << argv[1] << "' opened successfully.\n";
        while (std::getline(file, line)) {
            std::cout << line << '\n';
        }
    }
    return file;
}
