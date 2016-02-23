import abc
import os
import os.path
import xml.dom

import ArchiveComponent
import Info
import XmlUtils


class LabelMaker(XmlUtils.XmlUtils):
    """
    An abstract class of objects that can build PDS4 labels for
    ArchiveComponents.
    """

    def __init__(self, component, info):
        """
        Create the label for an ArchiveComponent with the help of a
        matching Info object to provide values for fields.
        """
        assert isinstance(component, ArchiveComponent.ArchiveComponent)
        self.component = component
        assert isinstance(info, Info.Info)
        self.info = info
        document = xml.dom.getDOMImplementation().createDocument(None,
                                                                 None,
                                                                 None)
        XmlUtils.XmlUtils.__init__(self, document)

        self.createDefaultXml()

    @abc.abstractmethod
    def createDefaultXml(self):
        """Create the XML label for the component."""
        pass

    @abc.abstractmethod
    def defaultXmlName(self):
        """The default name for the XML label for this type of component."""
        pass


def xmlSchemaCheck(filepath):
    """
    Test the XML label at the filepath against the PDS4 v1.5 XML
    schema, returning true iff it passes.
    """
    cmdTemplate = 'xmllint --noout --schema %s %s'
    exitCode = os.system(cmdTemplate %
                         ('./xml/PDS4_PDS_1500.xsd.xml', filepath))
    return exitCode == 0


def schematronCheck(filepath):
    """
    Test the XML label at the filepath against the PDS4 v1.5
    Schematron schema, returning true iff it passes.
    """
    cmdTemplate = 'java -jar probatron.jar %s %s' + \
        ' | xmllint -format -' + \
        ' | python Schematron.py'
    exitCode = os.system(cmdTemplate %
                         (filepath, './xml/PDS4_PDS_1500.sch.xml'))
    return exitCode == 0
