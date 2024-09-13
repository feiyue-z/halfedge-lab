#include <QCoreApplication>
#include <QCommandLineParser>
#include <QtCore>

#include <iostream>
#include <chrono>

#include "mesh.h"
#include "MeshDataStructure.h"

int main(int argc, char *argv[])
{
    QCoreApplication a(argc, argv);
    QCommandLineParser parser;
    parser.addHelpOption();
    parser.addPositionalArgument("config",  "Path of the config (.ini) file.");
    parser.process(a);

    // Check for invalid argument count
    const QStringList args = parser.positionalArguments();
    if (args.size() < 1) {
        std::cerr << "Not enough arguments. Please provide a path to a config file (.ini) as a command-line argument." << std::endl;
        a.exit(1);
        return 1;
    }

    // Parse common inputs
    QSettings settings( args[0], QSettings::IniFormat );
    QString infile  = settings.value("IO/infile").toString();
    QString outfile = settings.value("IO/outfile").toString();
    QString method  = settings.value("Method/method").toString();

    // A note about the representations of other parameters in the .ini files for the various methods:

    // args1:
    // Subdivide: number of iterations
    // Simplify:  number of faces to remove
    // Remesh:    number of iterations
    // Denoise:   number of iterations

    // args2:
    // Remesh: Tangential smoothing weight
    // Denoise: Smoothing parameter 1 (\Sigma_c)

    // args3:
    // Denoise: Smoothing parameter 2 (\Sigma_s)

    // args4:
    // Denoise: Kernel size (\rho)


    // Load
    Mesh m;
    m.loadFromFile(infile.toStdString());

    MeshDataStructure mesh;
    mesh.init(m);
    mesh.validate();

    // Start timing
    auto t0 = std::chrono::high_resolution_clock::now();

    // Switch on method
    if (method == "subdivide") {
        int numIterations = settings.value("Parameters/args1").toInt();

        mesh.subdivide(numIterations);
    } else if (method == "simplify") {
        int numFaceRemoves = settings.value("Parameters/args1").toInt();

        mesh.simplify(numFaceRemoves);
    } else if (method == "remesh") {
        int numIterations = settings.value("Parameters/args1").toInt();
        float weight = settings.value("Parameters/args2").toFloat();

        mesh.remesh(numIterations, weight);
    } else if (method == "noise") {
        int numIterations = settings.value("Parameters/args1").toInt();
        float scale = settings.value("Parameters/args2").toFloat();

        mesh.noise(numIterations,scale);
    } else if (method == "denoise") {
        int numIterations = settings.value("Parameters/args1").toInt();
        float sigma_c = settings.value("Parameters/args2").toFloat();
        float sigma_s = settings.value("Parameters/args3").toFloat();
        float rho = settings.value("Parameters/args4").toFloat();

        mesh.denoise(numIterations, sigma_c, sigma_s, rho);
    } else {
        std::cerr << "Error: Unknown method \"" << method.toUtf8().constData() << "\"" << std::endl;
    }

    // convert to .obj format
    m.initFromVectors(mesh.getVertices(), mesh.getFaces());

    // Finish timing
    auto t1 = std::chrono::high_resolution_clock::now();
    auto duration = duration_cast<std::chrono::milliseconds>(t1 - t0).count();
    std::cout << "Execution took " << duration << " milliseconds." << std::endl;

    // Save
    m.saveToFile(outfile.toStdString());

    a.exit();
}
