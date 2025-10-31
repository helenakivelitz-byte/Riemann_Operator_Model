# src/__init__.py
from .operator import RiemannOperator
OperatorBuilder = RiemannOperator  # back-compat alias

__all__ = ["RiemannOperator", "OperatorBuilder"]
