QT += widgets  printsupport widgets 
QT += core gui
requires(qtConfig(filedialog))



CONFIG += c++17 pthread #Quarter Coin



LIBS += -L/usr/local/lib64 -L/usr/local/lib




TARGET = symmetrizer
TEMPLATE = app


INCLUDEPATH += .
INCLUDEPATH    += /usr/local/includ
DEPENDPATH     += /usr/local/include




HEADERS       =  \
    matrix3x3.h \
    pointGroup.h \
    vector3.h



SOURCES       = main.cpp  \
    matrix3x3.cpp \
    pointGroup.cpp \
    vector3.cpp



DISTFILES +=








