"""Configuration management using pydantic-settings."""

from __future__ import annotations

from typing import Optional

from pydantic import model_validator
from pydantic_settings import BaseSettings


class Settings(BaseSettings):
    # Provider API keys (all optional — user may only use Ollama)
    openai_api_key: Optional[str] = None
    anthropic_api_key: Optional[str] = None

    # Ollama settings
    ollama_host: str = "http://localhost:11434"
    ollama_default_model: str = "qwen3.5:4b"
    ollama_keep_alive: str = "10m"

    # Default provider and model
    default_provider: str = "ollama"
    default_model: Optional[str] = None  # None means use provider default

    # LLM parameters
    default_temperature: float = 0.2
    default_max_tokens: int = 4096

    # Server settings
    server_host: str = "0.0.0.0"
    server_port: int = 8000

    # Code execution
    auto_execute_code: bool = True
    code_execution_timeout: int = 30

    # Memory compression
    memory_compression_threshold: int = 20  # Compress after this many messages
    memory_keep_recent: int = 4  # Keep last N messages intact when compressing

    # Logging
    log_level: str = "INFO"

    model_config = {"env_file": ".env", "env_file_encoding": "utf-8"}

    @model_validator(mode="before")
    @classmethod
    def _compat_and_clean_keys(cls, values: dict) -> dict:
        """Support plain OPENAI_API_KEY env var and strip surrounding quotes."""
        import os

        # Pick up plain OPENAI_API_KEY if openai_api_key not already set
        if not values.get("openai_api_key"):
            key = os.getenv("OPENAI_API_KEY", "")
            if key:
                values["openai_api_key"] = key

        # Strip surrounding quotes from API keys (common .env file artifact)
        for key_name in ("openai_api_key", "anthropic_api_key"):
            val = values.get(key_name)
            if isinstance(val, str):
                values[key_name] = val.strip("\"'")

        return values


_settings: Optional[Settings] = None


def get_settings() -> Settings:
    """Singleton accessor for application settings."""
    global _settings
    if _settings is None:
        _settings = Settings()
    return _settings
