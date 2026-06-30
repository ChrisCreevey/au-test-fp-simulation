#!/usr/bin/env python3
"""
Check if a set of taxa form a monophyletic clade in a tree.

A clade is monophyletic if:
  - All specified taxa are descendants of a single node
  - No non-specified taxa are descendants of that node
  - (i.e., the taxa form a complete, isolated subtree)
"""

import sys
import argparse
import re
from typing import Set, Tuple, Optional


class TreeNode:
    """Simple tree node class for Newick parsing."""

    def __init__(self, name: str = ""):
        self.name = name
        self.children = []
        self.parent = None

    def add_child(self, child):
        self.children.append(child)
        child.parent = self

    def is_leaf(self) -> bool:
        return len(self.children) == 0

    def get_all_descendants(self) -> Set[str]:
        """Get all leaf taxa under this node."""
        if self.is_leaf():
            return {self.name} if self.name else set()

        taxa = set()
        for child in self.children:
            taxa.update(child.get_all_descendants())
        return taxa


def parse_newick(newick_str: str) -> TreeNode:
    """
    Parse Newick format tree string into a tree structure.
    Simple parser that handles: (a,b)c; style trees.
    """
    newick_str = newick_str.strip().rstrip(';')

    stack = []
    current_tokens = []

    for char in newick_str:
        if char == '(':
            if current_tokens:
                # Save current token
                pass
            stack.append([])
            current_tokens = []
        elif char == ',':
            if current_tokens:
                token = ''.join(current_tokens).strip()
                if token and stack:
                    stack[-1].append(('token', token))
                current_tokens = []
        elif char == ')':
            if current_tokens:
                token = ''.join(current_tokens).strip()
                if token:
                    stack[-1].append(('token', token))
                current_tokens = []

            if stack:
                children_list = stack.pop()
                # Create internal node
                node = TreeNode()

                # Add children
                for item_type, item_value in children_list:
                    if item_type == 'token':
                        # Could be taxon or subtree (shouldn't happen in this simple parser)
                        child = TreeNode(item_value.split(':')[0])  # Remove branch length
                        node.add_child(child)
                    elif item_type == 'node':
                        node.add_child(item_value)

                # If this is the root
                if not stack:
                    return node
                else:
                    stack[-1].append(('node', node))
        else:
            current_tokens.append(char)

    # Handle any remaining tokens (shouldn't happen with valid Newick)
    return None


def parse_newick_simple(newick_str: str) -> TreeNode:
    """
    Parse Newick tree more robustly.
    Returns root node of the tree.
    """
    newick_str = newick_str.strip().rstrip(';')

    def parse_recursive(s: str, pos: int) -> Tuple[TreeNode, int]:
        """Parse from position pos, return (node, new_position)."""

        if s[pos] == '(':
            # Internal node
            node = TreeNode()
            pos += 1  # skip '('

            while s[pos] != ')':
                if s[pos] == ',':
                    pos += 1
                else:
                    child, pos = parse_recursive(s, pos)
                    node.add_child(child)

            pos += 1  # skip ')'

            # Parse optional label and branch length
            label = []
            while pos < len(s) and s[pos] not in '(),;':
                label.append(s[pos])
                pos += 1

            if label:
                label_str = ''.join(label)
                # Extract taxon name (before colon if present)
                if ':' in label_str:
                    node.name = label_str.split(':')[0]
                else:
                    node.name = label_str

            return node, pos
        else:
            # Leaf node
            label = []
            while pos < len(s) and s[pos] not in '(),;:':
                label.append(s[pos])
                pos += 1

            node = TreeNode(''.join(label))

            # Skip branch length
            if pos < len(s) and s[pos] == ':':
                pos += 1
                while pos < len(s) and s[pos] not in '(),;':
                    pos += 1

            return node, pos

    root, _ = parse_recursive(newick_str, 0)
    return root


def is_monophyletic(tree: TreeNode, target_taxa: Set[str]) -> bool:
    """
    Check if target_taxa form a monophyletic clade.

    Returns True if:
      - All target_taxa are present in the tree
      - They form a complete, isolated subtree
      - No non-target taxa are nested within them
    """

    if not target_taxa:
        return False

    # Get all taxa in tree
    all_taxa = tree.get_all_descendants()

    # Check if all target taxa exist
    if not target_taxa.issubset(all_taxa):
        return False

    # Find the MRCA (most recent common ancestor) of target taxa
    # by checking each node
    def find_mrca(node: TreeNode) -> Optional[TreeNode]:
        """Find the lowest node that has all target taxa as descendants."""
        node_taxa = node.get_all_descendants()

        if not target_taxa.issubset(node_taxa):
            return None

        # Check if this node is the MRCA
        # (all target taxa are here, and all descendants are target taxa)
        if node_taxa == target_taxa:
            return node

        # Check children
        mrcas = []
        for child in node.children:
            mrca = find_mrca(child)
            if mrca:
                mrcas.append(mrca)

        if mrcas:
            # Return the first one found (there should be only one)
            return mrcas[0]

        return None

    mrca = find_mrca(tree)

    if mrca is None:
        return False

    # Check if MRCA's descendants are exactly the target taxa
    mrca_taxa = mrca.get_all_descendants()

    return mrca_taxa == target_taxa


def main():
    parser = argparse.ArgumentParser(
        description='Check if specified taxa form a monophyletic clade in a tree'
    )
    parser.add_argument('--tree', required=True,
                        help='Newick format tree file')
    parser.add_argument('--taxa', required=True, nargs='+',
                        help='Taxon names to check for monophyly')

    args = parser.parse_args()

    # Read tree
    try:
        with open(args.tree, 'r') as f:
            tree_str = f.read().strip()
    except Exception as e:
        print(f"ERROR reading tree: {e}", file=sys.stderr)
        sys.exit(1)

    # Parse tree
    try:
        tree = parse_newick_simple(tree_str)
    except Exception as e:
        print(f"ERROR parsing tree: {e}", file=sys.stderr)
        sys.exit(1)

    if tree is None:
        print("ERROR: Failed to parse tree", file=sys.stderr)
        sys.exit(1)

    # Check monophyly
    target_taxa = set(args.taxa)
    is_mono = is_monophyletic(tree, target_taxa)

    # Output result as JSON for easy parsing
    import json
    result = {
        'monophyletic': is_mono,
        'taxa': sorted(target_taxa),
        'count': len(target_taxa)
    }

    print(json.dumps(result))
    sys.exit(0 if is_mono else 1)


if __name__ == '__main__':
    main()
