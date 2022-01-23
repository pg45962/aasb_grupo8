import unittest

from MyAlign import MyAlign
from MySeq import MySeq
from SubstMatrix import SubstMatrix
from AlignSeq import AlignSeq

class TestMySeq (unittest.TestCase):
    def transcricao (self):
        self.assertEqual(MySeq.transcricao('actg'),'ACUG')
        self.assertEqual(MySeq.transcricao('ACTG'),'ACUG')
    def compInverso (self):
        self.assertEqual(MySeq.compInverso('actg'),'CAGT')
        self.assertEqual(MySeq.compInverso('ACTG'),'CAGT')
    def traduzCodao (self):
        self.assertEqual(MySeq.traduzCodao('actg'),'T')
        self.assertEqual(MySeq.traduzCodao('ACTG'),'T')
    def validaER (self):
        self.assertEqual(MySeq.validaER('actg'),'Sequencia valida')
        self.assertEqual(MySeq.validaER('ACTG'),'Sequencia valida')
        self.assertEqual(MySeq.validaER('1234'),'Sequencia invalida')

class TestAlignSeq (unittest.TestCase):
    def smithWaterman (self):
        self.smithWaterman(AlignSeq.smithWaterman('WPCAWCWGPC','RFWP'),18)