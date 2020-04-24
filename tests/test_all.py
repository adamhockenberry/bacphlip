import pytest


def f():
    raise SystemExit(1)


def test_mytest():
    with pytest.raises(SystemExit):
        f()

def func(x):
    return x + 1


def test_answer():
    assert func(3) == 5


def test_imports():
    from Bio import SeqIO
    import bacphlip
