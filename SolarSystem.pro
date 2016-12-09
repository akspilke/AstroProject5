TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt


INCLUDEPATH += /usr/local/include
LIBS += -L/usr/local/lib -larmadillo

SOURCES += main.cpp \
    celestialbody.cpp \
    vec3.cpp \
    verlet.cpp \
    SolarSystem.cpp

HEADERS += \
    celestialbody.h \
    solarsystem.h \
    vec3.h \
    verlet.h

