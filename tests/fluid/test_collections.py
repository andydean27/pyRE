import pytest
from pyRE.fluid.models import Oil, Gas, Water
from pyRE.fluid.collections import FluidCollection

@pytest.fixture
def mock_fluids():
    """
    Fixture to provide mock fluid objects.
    """
    return {
        "oil": Oil(),  # Replace with appropriate initialization if needed
        "gas": Gas(),  # Replace with appropriate initialization if needed
        "water": Water()  # Replace with appropriate initialization if needed
    }

def test_initialization_with_valid_fluids(mock_fluids):
    """
    Test initialization with valid fluid models.
    """
    collection = FluidCollection(oil=mock_fluids["oil"], gas=mock_fluids["gas"], water=mock_fluids["water"])
    assert collection.oil == mock_fluids["oil"]
    assert collection.gas == mock_fluids["gas"]
    assert collection.water == mock_fluids["water"]

def test_initialization_with_invalid_fluid():
    """
    Test initialization with an invalid fluid type.
    """
    with pytest.raises(TypeError, match="Expected Oil model"):
        FluidCollection(oil="invalid_oil")  # Invalid type for oil

def test_add_fluid_valid(mock_fluids):
    """
    Test adding valid fluid models to the collection.
    """
    collection = FluidCollection()
    collection.add_fluid(mock_fluids["oil"])
    collection.add_fluid(mock_fluids["gas"])
    collection.add_fluid(mock_fluids["water"])

    assert collection.oil == mock_fluids["oil"]
    assert collection.gas == mock_fluids["gas"]
    assert collection.water == mock_fluids["water"]

def test_add_fluid_invalid():
    """
    Test adding an invalid fluid model to the collection.
    """
    collection = FluidCollection()
    with pytest.raises(TypeError, match="Expected Oil, Gas or Water model"):
        collection.add_fluid("invalid_fluid")  # Invalid type

def test_repr(mock_fluids):
    """
    Test the __repr__ method for proper string representation.
    """
    collection = FluidCollection(oil=mock_fluids["oil"], gas=mock_fluids["gas"], water=mock_fluids["water"])
    repr_str = repr(collection)
    assert "oil:" in repr_str
    assert "gas:" in repr_str
    assert "water:" in repr_str