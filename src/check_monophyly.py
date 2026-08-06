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


def is_monophyletic(tree: TreeNode, target_taxa: Set[str]) -> Optional[bool]:
    """
    Check if target_taxa form a "clan" - the unrooted analogue of a
    monophyletic clade (Wilkinson et al. 2007, "Of clades and clans").

    A rooted-descendant-set check (does some node's subtree equal exactly
    target_taxa?) is only correct if the tree happens to be rooted outside
    the target clade. ML tools like IQ-TREE write out arbitrarily rooted
    trees, so that assumption can't be relied on: if the root falls inside
    the true clade, no single node's descendant set will equal target_taxa
    even though the taxa are still separated from everything else by one
    edge in the underlying unrooted tree.

    Instead, we test every internal edge's bipartition (a node's descendant
    set vs. its complement). Both sides of a bipartition are checked because
    rerooting can only reshuffle the side of the tree containing the root -
    the side that doesn't contain the root is always left as a clean
    subtree, so it always shows up as either a node's descendant set or its
    complement. This makes the test rooting-invariant.

    Returns:
      - True if the (present) target taxa form a clan
      - False if they are violated (found in the tree but not grouped together)
      - None if fewer than 2 target taxa are present in the tree (untestable,
        mirroring the "?" result in clan_check)
    """

    if not target_taxa:
        return False

    all_taxa = tree.get_all_descendants()

    # Only test the subset of target taxa actually present in this tree
    present_target = target_taxa & all_taxa
    if len(present_target) < 2:
        return None

    def matches_split(node: TreeNode) -> bool:
        node_taxa = node.get_all_descendants()
        if node_taxa == present_target:
            return True
        if (all_taxa - node_taxa) == present_target:
            return True
        return False

    def walk(node: TreeNode) -> bool:
        if matches_split(node):
            return True
        return any(walk(child) for child in node.children)

    return walk(tree)


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
    # Exit codes: 0 = clan present, 1 = clan violated, 2 = untestable (< 2 taxa present)
    sys.exit(0 if is_mono is True else (2 if is_mono is None else 1))


if __name__ == '__main__':
    main()
