import pytest
from molbench.tree import Node, DummyNode


def test_node_construction_links_parent():
    root = Node("root")
    child = Node("child", root)
    assert child in root.children
    assert child.parent is root


def test_node_no_parent_is_root():
    n = Node("x")
    assert n.is_root
    assert n.parent is None


def test_node_with_children_not_leaf():
    root = Node("root")
    Node("child", root)
    assert not root.is_leave


def test_node_no_children_is_leaf():
    n = Node("leaf")
    assert n.is_leave


def test_traverse_depth_first_order():
    root = Node("r")
    a = Node("a", root)
    b = Node("b", root)
    Node("a1", a)
    Node("a2", a)
    Node("b1", b)
    values = [n.value for n in root.traverse()]
    # DFS: r, a, a1, a2, b, b1
    assert values == ["r", "a", "a1", "a2", "b", "b1"]


def test_traverse_breadth_first_skips_dummy():
    root = DummyNode()
    a = Node("a", root)
    Node("a1", a)
    Node("b", root)
    values = [n.value for n in root.traverse_breadth_first()]
    assert values == ["a", "b", "a1"]


def test_traverse_generations_two_levels():
    root = Node("root")
    a = Node("a", root)
    Node("b", root)
    Node("a1", a)
    gens = list(root.traverse_generations())
    assert len(gens) == 3
    assert tuple(n.value for n in gens[0]) == ("root",)
    assert set(n.value for n in gens[1]) == {"a", "b"}
    assert tuple(n.value for n in gens[2]) == ("a1",)


def test_traverse_generations_skips_dummy():
    root = DummyNode()
    a = Node("a", root)
    Node("a1", a)
    gens = list(root.traverse_generations())
    # DummyNode not yielded as first generation
    assert len(gens) == 2
    assert tuple(n.value for n in gens[0]) == ("a",)


def test_walk_leaves_only_terminal():
    root = Node("root")
    a = Node("a", root)
    b = Node("b", root)
    leaf_a = Node("a1", a)
    leaf_b = Node("b1", b)
    leaves = list(root.walk_leaves())
    assert set(n.value for n in leaves) == {"a1", "b1"}
    assert a not in leaves
    assert root not in leaves


def test_sort_recursive():
    root = Node("root")
    Node("c", root)
    Node("a", root)
    Node("b", root)
    root.sort(key=lambda n: n.value)
    assert [n.value for n in root.children] == ["a", "b", "c"]


def test_path_to_root():
    root = Node("root")
    mid = Node("mid", root)
    leaf = Node("leaf", mid)
    path = leaf.path_to_root()
    assert [n.value for n in path] == ["leaf", "mid", "root"]


def test_path_to_root_excludes_dummy():
    dummy = DummyNode()
    child = Node("child", dummy)
    leaf = Node("leaf", child)
    path = leaf.path_to_root()
    assert "DummyNode" not in [str(n) for n in path]
    assert [n.value for n in path] == ["leaf", "child"]


def test_depth_leaf():
    n = Node("x")
    assert n.depth() == 0


def test_depth_chain():
    root = Node("root")
    mid = Node("mid", root)
    Node("leaf", mid)
    assert root.depth() == 2
    assert mid.depth() == 1


def test_width_leaf():
    n = Node("x")
    assert n.width() == 1


def test_width_branching():
    root = Node("root")
    Node("a", root)
    Node("b", root)
    # root.width = a.width + b.width + 1 = 1 + 1 + 1 = 3
    assert root.width() == 3


def test_dummy_node_traverse_skips_self():
    dummy = DummyNode()
    Node("a", dummy)
    Node("b", dummy)
    values = [n.value for n in dummy.traverse()]
    assert values == ["a", "b"]


def test_to_string_default():
    n = Node("key")
    assert n.to_string(42) == "42"


def test_to_string_custom():
    n = Node("key", to_string=lambda x: f"val={x}")
    assert n.to_string(5) == "val=5"


def test_node_str():
    n = Node("hello")
    assert str(n) == "Node(hello)"


def test_dummy_node_str():
    d = DummyNode()
    assert str(d) == "DummyNode"
