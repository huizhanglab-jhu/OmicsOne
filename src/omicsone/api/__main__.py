def main():
    import uvicorn

    uvicorn.run(
        "omicsone.api.app:app",
        host="127.0.0.1",
        port=8001,
        reload=False,
    )


if __name__ == "__main__":
    main()
