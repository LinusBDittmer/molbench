from collections import deque
from itertools import chain


class Node:
    def __init__(self, value, parent: 'Node' = None,
                 to_string: callable = None) -> None:
        self.value = value
        self.children: list['Node'] = []
        self.parent: 'Node' | None = parent
        if parent is not None:
            parent.children.append(self)
        # callable to convert the found keys to string
        self.to_string = lambda x: str(x) if to_string is None else to_string

    def __str__(self):
        return f"Node({self.value})"

    @property
    def is_leave(self) -> bool:
        return not self.children

    @property
    def is_root(self) -> bool:
        return self.parent is None

    def traverse(self):
        # depth first algorithm
        yield self
        for child in self.children:
            for n in child.traverse():
                yield n

    def traverse_breadth_first(self):
        # skip all DummyNodes
        queue = deque()
        if isinstance(self, DummyNode):
            queue.extend(self.children)
        else:
            queue.append(self)
        while len(queue) > 0:
            node = queue.popleft()
            if not isinstance(node, DummyNode):
                yield node
            queue.extend(node.children)

    def traverse_generations(self):
        nodes = (self,)
        if not isinstance(self, DummyNode):
            yield nodes
        while any(node.children for node in nodes):
            nodes = tuple(chain.from_iterable(
                (n for n in node.children if not isinstance(n, DummyNode))
                for node in nodes
            ))
            if nodes:  # there are non dummy nodes in the generation
                yield nodes

    def walk_leaves(self):
        for n in self.traverse():
            if not n.children:
                yield n

    def sort(self, **kwargs):
        # sort the tree by value
        self.children.sort(**kwargs)
        for child in self.children:
            child.sort(**kwargs)

    def path_to_root(self) -> list['Node']:
        # path to the root node including self and excluding DummyNodes
        path: list[Node] = []
        if not isinstance(self, DummyNode):
            path.append(self)
        node = self
        while not node.is_root:
            node = node.parent
            if not isinstance(node, DummyNode):
                path.append(node)
        return path

    def depth(self) -> int:
        if self.children:
            return max(n.depth() for n in self.children) + 1
        else:
            return 0

    def width(self) -> int:
        if self.children:
            return sum(n.width() for n in self.children) + 1
        else:
            return 1


class DummyNode(Node):
    def __init__(self, parent: 'Node' = None) -> None:
        self.children = []
        self.parent: 'Node' | None = parent

    def __str__(self):
        return "DummyNode"

    def traverse(self):
        # skip dummy nodes
        for child in self.children:
            for n in child.traverse():
                yield n
