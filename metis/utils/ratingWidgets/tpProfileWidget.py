from PySide6 import QtWidgets


class TPProfile(QtWidgets.QWidget):
    def __init__(self, settings):
        super(TPProfile, self).__init__()

        # Layout
        self.mainlayout = QtWidgets.QVBoxLayout()

        # Widgets
        introductionText = QtWidgets.QLabel()
        introductionText.setText(settings["introText"])
        introductionText.setWordWrap(True)

        propertyList = QtWidgets.QLabel()
        propertyList.setText(
            "<html><ul>"
            + "".join([f"<li>{name}</li>" for name in settings["propertyLabels"]])
            + "</ul></html>"
        )

        self.displayTPP = QtWidgets.QLabel()
        self.displayTPP.setText(
            "<html><ul>"
            + "".join([f"<li>{name}</li>" for name in settings["propertyLabels"]])
            + "</ul></html>"
        )

        self.mainlayout.addWidget(introductionText)
        self.mainlayout.addWidget(propertyList)
        self.mainlayout.addWidget(QHLine())
        self.mainlayout.addWidget(self.displayTPP)
        self.setLayout(self.mainlayout)

    def updateMolProperties(self, propertyDict):
        "update displayed properties when molecule is changed"
        stringToDisplay = ["<html><table>"]
        for i in propertyDict:
            stringToDisplay += f"<tr><td>{i}:</td> <td align='right'>{(propertyDict[i]*100):.2f}%</td> </tr>"
        stringToDisplay += "</table></html>"
        self.displayTPP.setText("".join(stringToDisplay))


class QHLine(QtWidgets.QFrame):
    def __init__(self):
        super(QHLine, self).__init__()
        self.setFrameShape(QtWidgets.QFrame.HLine)
        self.setFrameShadow(QtWidgets.QFrame.Sunken)


class QVLine(QtWidgets.QFrame):
    def __init__(self):
        super(QVLine, self).__init__()
        self.setFrameShape(QtWidgets.QFrame.VLine)
        self.setFrameShadow(QtWidgets.QFrame.Sunken)
